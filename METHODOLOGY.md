# Promoter Coordinate Correction Methodology

## Overview

This tool addresses the problem of inconsistent promoter region coordinates in BED files by realigning them with authoritative gene coordinates (gencode/refseq) using transcription start site (TSS) data from CAGE experiments and chromatin marks from chromHMM/histone modifications.

## Problem Statement

Promoter regions defined in BED files often suffer from:
- Inconsistent coordinate systems between different databases
- Outdated gene annotations
- Arbitrary upstream/downstream boundaries
- Lack of experimental validation

## Methodology

### 1. Data Integration

The tool integrates multiple genomic datasets:

- **Original Promoter BED**: The inconsistent promoter coordinates to be corrected
- **Gene Coordinates**: Authoritative gene boundaries from gencode or refseq (GTF/BED format)
- **CAGE TSS Data**: Experimentally determined transcription start sites from FANTOM project
- **Histone Marks**: Chromatin modifications indicating active promoters (H3K4me3, H3K27ac)

### 2. Coordinate Correction Algorithm

#### Step 1: TSS-Based Realignment
For each promoter region:
1. Find the nearest experimentally validated TSS from CAGE data
2. Calculate new promoter boundaries relative to the TSS:
   - **Positive strand**: TSS - upstream_bp to TSS + downstream_bp
   - **Negative strand**: TSS - downstream_bp to TSS + upstream_bp
   - Default: 1000bp upstream, 500bp downstream

#### Step 2: Histone Mark Validation
1. Check overlap between corrected coordinates and histone marks
2. Validate promoter activity using:
   - **H3K4me3**: Active promoter mark
   - **H3K27ac**: Active enhancer/promoter mark
3. Require minimum overlap threshold (default: 30%)

#### Step 3: Quality Control
1. Compare original vs. corrected coordinates
2. Flag significant changes for manual review
3. Maintain regions without nearby TSS with original coordinates
4. Report statistics on correction success rates

### 3. Output Format

The corrected BED file contains:
- **chr**: Chromosome name
- **start**: Corrected start coordinate
- **end**: Corrected end coordinate  
- **name**: Promoter identifier (with "_corrected" suffix)
- **score**: TSS confidence score from CAGE data
- **strand**: Gene orientation

## Usage

### Basic Usage
```bash
python promoter_corrector.py \
    -p original_promoters.bed \
    -g gencode_genes.gtf \
    -t cage_tss.bed \
    --histone H3K4me3_marks.bed \
    -o corrected_promoters.bed
```

### Advanced Options
```bash
python promoter_corrector.py \
    -p promoters.bed \
    -g refseq_genes.bed \
    --gene-source refseq \
    -t fantom_tss.bed \
    --histone H3K4me3.bed \
    --histone-type H3K4me3 \
    --upstream 2000 \
    --downstream 500 \
    --min-overlap 0.4 \
    -o corrected_promoters.bed \
    --verbose
```

## File Format Requirements

### Input Files

#### 1. Original Promoter BED
Standard BED format with at least 3 columns:
```
chr1    1000    2500    promoter_1    0    +
chr1    5000    6500    promoter_2    0    -
```

#### 2. Gene Coordinates
**GTF format** (preferred for gencode):
```
chr1    HAVANA    gene    1000    5000    .    +    .    gene_id "ENSG00000001"; gene_name "GENE1";
```

**BED format** (alternative):
```
chr1    1000    5000    ENSG00000001    0    +
```

#### 3. CAGE TSS Data
BED format with TSS positions:
```
chr1    1500    1501    tss_1    100    +
chr1    3000    3001    tss_2    150    -
```

#### 4. Histone Mark Data
BED format with chromatin modification regions:
```
chr1    1200    1800    H3K4me3_peak_1    500    .
chr1    2900    3200    H3K4me3_peak_2    750    .
```

### Output File

Corrected promoter BED with metadata:
```
# Corrected promoter coordinates
# Columns: chr, start, end, name, score, strand
chr1    500    2000    promoter_1_corrected    100    +
chr1    2500    4000    promoter_2_corrected    150    -
```

## Validation Strategy

### 1. Distance Metrics
- Calculate distance between original and corrected coordinates
- Flag promoters with >1kb shifts for manual review

### 2. Overlap Analysis
- Measure histone mark coverage in corrected regions
- Compare with original regions to assess improvement

### 3. Gene Association
- Verify promoter-gene associations remain biologically relevant
- Check for promoter-TSS distance consistency

## Quality Metrics

The tool reports:
- **Total promoters processed**
- **Successfully corrected** (with nearby TSS)
- **Histone-validated** (with sufficient chromatin marks)
- **No TSS found** (kept original coordinates)
- **No correction needed** (minimal coordinate change)

## Biological Rationale

### TSS-Centric Definition
- Promoters are regulatory regions surrounding transcription start sites
- CAGE data provides experimental evidence for active TSS positions
- More accurate than computational predictions alone

### Chromatin State Validation
- H3K4me3: Marks active and poised promoters
- H3K27ac: Indicates active regulatory elements
- Provides functional validation of predicted promoter boundaries

### Strand-Aware Boundaries
- Accounts for gene orientation in promoter definition
- Ensures biological relevance of upstream/downstream regions

## Limitations and Considerations

1. **Cell Type Specificity**: CAGE and histone data may be cell-type specific
2. **Alternative Promoters**: Genes may have multiple promoters not captured
3. **Enhancer Confusion**: Strong enhancers may be misidentified as promoters
4. **Assembly Differences**: Coordinate systems must match between datasets

## Best Practices

1. **Data Consistency**: Ensure all input files use the same genome assembly
2. **Parameter Tuning**: Adjust upstream/downstream distances based on organism
3. **Manual Review**: Examine flagged promoters with large coordinate changes
4. **Validation**: Cross-reference with independent datasets when available

## Example Workflow

1. Download CAGE TSS data from FANTOM5
2. Obtain gene coordinates from gencode/refseq
3. Download chromHMM states or histone ChIP-seq peaks
4. Run promoter correction tool
5. Validate results using expression data or additional chromatin marks
6. Apply corrected coordinates to downstream analyses

This methodology provides a systematic approach to improve promoter coordinate accuracy using multiple lines of experimental evidence.