# Promoter Coordinate Correction Tool

[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A bioinformatics tool to correct inconsistent promoter BED file coordinates by aligning them with authoritative gene coordinates using CAGE TSS data and chromatin modifications.

## 🧬 Problem

Promoter regions in genomic datasets often have inconsistent coordinates due to:
- Different annotation versions
- Arbitrary boundary definitions  
- Lack of experimental validation
- Inconsistent coordinate systems

## 🔬 Solution

This tool provides a systematic methodology to correct promoter coordinates by:

1. **TSS-based realignment** using CAGE data from FANTOM
2. **Chromatin state validation** using histone marks (H3K4me3, H3K27ac)
3. **Gene boundary integration** from gencode/refseq annotations
4. **Quality control metrics** for validation

## 📁 Files

- `promoter_corrector.py` - Main correction tool
- `METHODOLOGY.md` - Detailed methodology and biological rationale
- `test_promoter_corrector.py` - Test script with sample data
- `requirements.txt` - Python dependencies

## ⚡ Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/arnabmukho/arnabmukho.git
cd arnabmukho

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```bash
python promoter_corrector.py \
    -p original_promoters.bed \
    -g gencode_genes.gtf \
    -t cage_tss.bed \
    --histone H3K4me3_marks.bed \
    -o corrected_promoters.bed
```

### Test the Tool

```bash
python test_promoter_corrector.py
```

## 📋 Input File Formats

### Original Promoter BED
```
chr1    1000    2500    promoter_1    0    +
chr1    5000    6500    promoter_2    0    -
```

### Gene Coordinates (GTF/BED)
```
chr1    HAVANA    gene    1000    5000    .    +    .    gene_id "ENSG00000001";
```

### CAGE TSS Data
```
chr1    1500    1501    tss_1    100    +
chr1    3000    3001    tss_2    150    -
```

### Histone Marks
```
chr1    1200    1800    H3K4me3_peak_1    500    .
chr1    2900    3200    H3K4me3_peak_2    750    .
```

## 🎛️ Command Line Options

```bash
python promoter_corrector.py [OPTIONS]

Required:
  -p, --promoters     Original promoter BED file
  -g, --genes         Gene coordinates (GTF/BED)
  -t, --tss          CAGE TSS data (BED)
  -o, --output       Output corrected BED file

Optional:
  --histone          Histone mark data (BED)
  --gene-source      gencode|refseq (default: gencode)
  --histone-type     H3K4me3|H3K27ac (default: H3K4me3)
  --upstream         Upstream bp from TSS (default: 1000)
  --downstream       Downstream bp from TSS (default: 500)
  --min-overlap      Min histone overlap fraction (default: 0.3)
  --verbose          Enable verbose logging
```

## 🔬 Methodology Overview

### 1. TSS-Based Correction
- Find nearest experimentally validated TSS for each promoter
- Recalculate boundaries relative to TSS position
- Account for gene strand orientation

### 2. Chromatin Validation
- Validate corrected regions using histone marks
- Require minimum overlap with H3K4me3/H3K27ac
- Flag regions without chromatin support

### 3. Quality Control
- Report correction statistics
- Identify promoters requiring manual review
- Maintain original coordinates when no TSS found

## 📊 Output Statistics

The tool reports:
- **Total promoters processed**
- **Successfully corrected** (with nearby TSS)
- **Histone-validated** (with chromatin marks)
- **No TSS found** (original coordinates kept)
- **No correction needed** (minimal change)

## 🧪 Example Output

```
2024-01-01 10:00:00 - INFO - Loaded 1000 original promoter regions
2024-01-01 10:00:01 - INFO - Loaded 20000 genes from gencode
2024-01-01 10:00:02 - INFO - Loaded 15000 TSS regions from CAGE data
2024-01-01 10:00:03 - INFO - Loaded 8000 H3K4me3 regions
2024-01-01 10:00:05 - INFO - Promoter correction completed!
2024-01-01 10:00:05 - INFO - Total promoters: 1000
2024-01-01 10:00:05 - INFO - TSS-corrected: 850
2024-01-01 10:00:05 - INFO - Histone-validated: 720
2024-01-01 10:00:05 - INFO - No TSS found: 100
2024-01-01 10:00:05 - INFO - No correction needed: 50
```

## 📚 Data Sources

### Recommended Datasets

1. **CAGE TSS Data**: [FANTOM5](https://fantom.gsc.riken.jp/)
2. **Gene Coordinates**: [GENCODE](https://www.gencodegenes.org/) or [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)
3. **Histone Marks**: [ENCODE](https://www.encodeproject.org/) or [Roadmap Epigenomics](http://www.roadmapepigenomics.org/)
4. **chromHMM States**: [Roadmap chromHMM](http://www.roadmapepigenomics.org/)

### Data Processing Tips

1. Ensure all datasets use the same genome assembly (e.g., hg38)
2. Filter CAGE peaks by expression threshold
3. Use cell-type specific histone marks when available
4. Consider multiple histone modifications for validation

## 🔧 Advanced Usage

### Custom Promoter Boundaries
```bash
python promoter_corrector.py \
    --upstream 2000 \
    --downstream 500 \
    -p promoters.bed -g genes.gtf -t tss.bed -o output.bed
```

### Multiple Histone Mark Validation
```bash
# Run with H3K4me3
python promoter_corrector.py --histone H3K4me3.bed --histone-type H3K4me3 [other options]

# Run with H3K27ac
python promoter_corrector.py --histone H3K27ac.bed --histone-type H3K27ac [other options]
```

### RefSeq Gene Coordinates
```bash
python promoter_corrector.py \
    --gene-source refseq \
    -g refseq_genes.bed [other options]
```

## 🐛 Troubleshooting

### Common Issues

1. **No TSS found**: Increase search distance or check chromosome naming consistency
2. **Low histone validation**: Adjust `--min-overlap` threshold or check mark quality
3. **Coordinate mismatches**: Verify all files use same genome assembly
4. **Memory issues**: Process chromosomes separately for large datasets

### File Format Issues

- Ensure BED files are tab-separated
- Check chromosome naming convention (chr1 vs 1)
- Verify coordinate systems (0-based vs 1-based)
- Remove header lines from data files

## 📄 Citation

If you use this tool in your research, please cite:

```
Promoter Coordinate Correction Tool
Author: Arnab Mukhopadhyay
GitHub: https://github.com/arnabmukho/arnabmukho
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📞 Contact

- **Author**: Arnab Mukhopadhyay
- **Email**: arnabbiotech.gen@gmail.com
- **LinkedIn**: [arnabmukho](https://www.linkedin.com/in/arnabmukho)

---

<p align="center">
  <img src="https://readme-typing-svg.demolab.com?font=Fira+Code&pause=2000&color=7ED957&width=435&lines=Fixing+promoter+coordinates%2C+one+TSS+at+a+time!+%F0%9F%A7%AC" alt="Animated typing" />
</p>