#!/usr/bin/env python3
"""
Test script for the promoter corrector tool.
Creates sample data and tests the functionality.
"""

import os
import tempfile
import sys
from promoter_corrector import PromoterCorrector

def create_sample_data():
    """Create sample data files for testing."""
    
    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix="promoter_test_")
    
    # Sample promoter BED file (inconsistent coordinates)
    promoter_file = os.path.join(temp_dir, "original_promoters.bed")
    with open(promoter_file, 'w') as f:
        f.write("chr1\t1000\t2500\tpromoter_1\t0\t+\n")
        f.write("chr1\t5000\t6500\tpromoter_2\t0\t-\n")
        f.write("chr2\t10000\t11500\tpromoter_3\t0\t+\n")
    
    # Sample gene coordinates (GTF-like format converted to BED)
    gene_file = os.path.join(temp_dir, "genes.bed")
    with open(gene_file, 'w') as f:
        f.write("chr1\t1500\t4000\tENSG00000001\t0\t+\n")
        f.write("chr1\t5500\t8000\tENSG00000002\t0\t-\n")
        f.write("chr2\t9500\t12000\tENSG00000003\t0\t+\n")
    
    # Sample CAGE TSS data
    tss_file = os.path.join(temp_dir, "cage_tss.bed")
    with open(tss_file, 'w') as f:
        f.write("chr1\t1600\t1601\ttss_1\t100\t+\n")
        f.write("chr1\t5600\t5601\ttss_2\t150\t-\n")
        f.write("chr2\t10100\t10101\ttss_3\t120\t+\n")
    
    # Sample histone mark data (H3K4me3)
    histone_file = os.path.join(temp_dir, "H3K4me3.bed")
    with open(histone_file, 'w') as f:
        f.write("chr1\t1400\t1800\tH3K4me3_peak_1\t500\t.\n")
        f.write("chr1\t5400\t5800\tH3K4me3_peak_2\t750\t.\n")
        f.write("chr2\t9900\t10300\tH3K4me3_peak_3\t600\t.\n")
    
    return temp_dir, promoter_file, gene_file, tss_file, histone_file

def test_promoter_correction():
    """Test the promoter correction functionality."""
    
    print("Creating sample data...")
    temp_dir, promoter_file, gene_file, tss_file, histone_file = create_sample_data()
    
    print("Initializing promoter corrector...")
    corrector = PromoterCorrector(upstream_bp=1000, downstream_bp=500)
    
    print("Loading data files...")
    # Load original promoters
    corrector.read_original_promoters(promoter_file)
    
    # Load gene coordinates  
    corrector.read_gene_coordinates(gene_file, "test")
    
    # Load TSS data
    corrector.read_cage_tss_data(tss_file)
    
    # Load histone marks
    corrector.read_histone_marks(histone_file, "H3K4me3")
    
    print("Performing coordinate correction...")
    output_file = os.path.join(temp_dir, "corrected_promoters.bed")
    corrected_promoters, stats = corrector.correct_promoter_coordinates(output_file)
    
    print("\nCorrection Results:")
    print(f"Total promoters: {stats['total']}")
    print(f"TSS-corrected: {stats['tss_corrected']}")
    print(f"Histone-validated: {stats['histone_validated']}")
    print(f"No TSS found: {stats['no_tss_found']}")
    print(f"No correction needed: {stats['no_correction_needed']}")
    
    print(f"\nOutput written to: {output_file}")
    
    # Show corrected coordinates
    print("\nCorrected Promoter Coordinates:")
    for promoter in corrected_promoters:
        print(f"{promoter['name']}: {promoter['chr']}:{promoter['start']}-{promoter['end']}")
        if 'original_start' in promoter:
            print(f"  Original: {promoter['chr']}:{promoter['original_start']}-{promoter['original_end']}")
            print(f"  TSS distance: {promoter.get('tss_distance', 'N/A')}")
            print(f"  Histone validated: {promoter.get('histone_validated', False)}")
    
    print(f"\nTest completed successfully!")
    print(f"Temporary files in: {temp_dir}")
    
    return True

if __name__ == "__main__":
    try:
        test_promoter_correction()
    except Exception as e:
        print(f"Test failed with error: {e}")
        sys.exit(1)