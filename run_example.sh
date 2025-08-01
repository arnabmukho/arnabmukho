#!/bin/bash
# Example script to demonstrate promoter coordinate correction

echo "Promoter Coordinate Correction Example"
echo "====================================="

# Create example data directory
mkdir -p example_data

# Create sample promoter BED file
echo "Creating sample promoter BED file..."
cat > example_data/promoters.bed << EOF
chr1	1000	2500	promoter_1	0	+
chr1	5000	6500	promoter_2	0	-
chr2	10000	11500	promoter_3	0	+
chr3	15000	16500	promoter_4	0	+
EOF

# Create sample gene coordinates
echo "Creating sample gene coordinates..."
cat > example_data/genes.bed << EOF
chr1	1500	4000	ENSG00000001	0	+
chr1	5500	8000	ENSG00000002	0	-
chr2	9500	12000	ENSG00000003	0	+
chr3	14500	17000	ENSG00000004	0	+
EOF

# Create sample CAGE TSS data
echo "Creating sample CAGE TSS data..."
cat > example_data/cage_tss.bed << EOF
chr1	1600	1601	tss_1	100	+
chr1	5600	5601	tss_2	150	-
chr2	10100	10101	tss_3	120	+
chr3	15100	15101	tss_4	130	+
EOF

# Create sample histone mark data
echo "Creating sample histone mark data..."
cat > example_data/H3K4me3.bed << EOF
chr1	1400	1800	H3K4me3_peak_1	500	.
chr1	5400	5800	H3K4me3_peak_2	750	.
chr2	9900	10300	H3K4me3_peak_3	600	.
chr3	14900	15300	H3K4me3_peak_4	650	.
EOF

echo "Running promoter coordinate correction..."
python promoter_corrector.py \
    -p example_data/promoters.bed \
    -g example_data/genes.bed \
    -t example_data/cage_tss.bed \
    --histone example_data/H3K4me3.bed \
    -o example_data/corrected_promoters.bed \
    --verbose

echo ""
echo "Correction completed! Results saved to example_data/corrected_promoters.bed"
echo ""
echo "Original vs Corrected coordinates:"
echo "=================================="

echo "Original promoters:"
cat example_data/promoters.bed

echo ""
echo "Corrected promoters:"
cat example_data/corrected_promoters.bed

echo ""
echo "Example completed successfully!"