#!/usr/bin/env python3
"""
Promoter Coordinate Correction Tool

This script corrects inconsistent promoter BED file coordinates by aligning them 
with gencode/refseq gene coordinates using chromHMM/histone marks and CAGE TSS data.

Author: Arnab Mukhopadhyay
"""

import argparse
import pandas as pd
import numpy as np
import sys
from collections import defaultdict
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class PromoterCorrector:
    """
    A class to correct promoter coordinates using TSS data and histone marks.
    """
    
    def __init__(self, upstream_bp=1000, downstream_bp=500):
        """
        Initialize the promoter corrector.
        
        Args:
            upstream_bp (int): Base pairs upstream of TSS to include in promoter (default: 1000)
            downstream_bp (int): Base pairs downstream of TSS to include in promoter (default: 500)
        """
        self.upstream_bp = upstream_bp
        self.downstream_bp = downstream_bp
        self.genes = {}
        self.tss_sites = {}
        self.histone_marks = {}
        self.original_promoters = {}
        
    def read_bed_file(self, filename, file_type="bed"):
        """
        Read a BED format file.
        
        Args:
            filename (str): Path to BED file
            file_type (str): Type of file being read
            
        Returns:
            pandas.DataFrame: BED file contents
        """
        try:
            # Standard BED file has at least 3 columns: chr, start, end
            df = pd.read_csv(filename, sep='\t', header=None, 
                           names=['chr', 'start', 'end', 'name', 'score', 'strand'],
                           dtype={'chr': str, 'start': int, 'end': int})
            logger.info(f"Read {len(df)} entries from {file_type} file: {filename}")
            return df
        except Exception as e:
            logger.error(f"Error reading {file_type} file {filename}: {e}")
            return None
            
    def read_gene_coordinates(self, filename, source="gencode"):
        """
        Read gene coordinates from gencode or refseq file.
        
        Args:
            filename (str): Path to gene coordinate file (GTF or BED format)
            source (str): Source of gene data ("gencode" or "refseq")
        """
        try:
            if filename.endswith('.gtf') or filename.endswith('.gff'):
                # GTF format
                df = pd.read_csv(filename, sep='\t', header=None, comment='#',
                               names=['chr', 'source', 'feature', 'start', 'end', 
                                     'score', 'strand', 'frame', 'attributes'])
                # Filter for gene entries
                genes_df = df[df['feature'] == 'gene'].copy()
                
                # Extract gene information from attributes
                for idx, row in genes_df.iterrows():
                    attrs = dict(item.strip().split(' ', 1) for item in row['attributes'].split(';') if item.strip())
                    gene_id = attrs.get('gene_id', '').strip('"')
                    gene_name = attrs.get('gene_name', '').strip('"')
                    
                    if gene_id:
                        self.genes[gene_id] = {
                            'chr': row['chr'],
                            'start': row['start'],
                            'end': row['end'],
                            'strand': row['strand'],
                            'name': gene_name or gene_id,
                            'source': source
                        }
            else:
                # BED format
                df = self.read_bed_file(filename, f"{source} genes")
                if df is not None:
                    for idx, row in df.iterrows():
                        gene_id = row.get('name', f"gene_{idx}")
                        self.genes[gene_id] = {
                            'chr': row['chr'],
                            'start': row['start'],
                            'end': row['end'],
                            'strand': row.get('strand', '+'),
                            'name': gene_id,
                            'source': source
                        }
                        
            logger.info(f"Loaded {len(self.genes)} genes from {source}")
            
        except Exception as e:
            logger.error(f"Error reading gene coordinates from {filename}: {e}")
            
    def read_cage_tss_data(self, filename):
        """
        Read CAGE TSS data from FANTOM.
        
        Args:
            filename (str): Path to CAGE TSS file (BED format)
        """
        try:
            df = self.read_bed_file(filename, "CAGE TSS")
            if df is not None:
                for idx, row in df.iterrows():
                    tss_id = row.get('name', f"tss_{idx}")
                    chr_tss = f"{row['chr']}:{row['start']}-{row['end']}"
                    
                    if chr_tss not in self.tss_sites:
                        self.tss_sites[chr_tss] = []
                        
                    self.tss_sites[chr_tss].append({
                        'chr': row['chr'],
                        'start': row['start'],
                        'end': row['end'],
                        'strand': row.get('strand', '+'),
                        'name': tss_id,
                        'score': row.get('score', 0)
                    })
                    
            logger.info(f"Loaded {len(self.tss_sites)} TSS regions from CAGE data")
            
        except Exception as e:
            logger.error(f"Error reading CAGE TSS data from {filename}: {e}")
            
    def read_histone_marks(self, filename, mark_type="H3K4me3"):
        """
        Read chromHMM or histone mark data.
        
        Args:
            filename (str): Path to histone mark file (BED format)
            mark_type (str): Type of histone mark (e.g., "H3K4me3", "H3K27ac")
        """
        try:
            df = self.read_bed_file(filename, f"histone mark {mark_type}")
            if df is not None:
                mark_regions = []
                for idx, row in df.iterrows():
                    mark_regions.append({
                        'chr': row['chr'],
                        'start': row['start'],
                        'end': row['end'],
                        'score': row.get('score', 0)
                    })
                    
                self.histone_marks[mark_type] = mark_regions
                logger.info(f"Loaded {len(mark_regions)} {mark_type} regions")
                
        except Exception as e:
            logger.error(f"Error reading histone mark data from {filename}: {e}")
            
    def read_original_promoters(self, filename):
        """
        Read the original promoter BED file that needs correction.
        
        Args:
            filename (str): Path to original promoter BED file
        """
        try:
            df = self.read_bed_file(filename, "original promoters")
            if df is not None:
                for idx, row in df.iterrows():
                    promoter_id = row.get('name', f"promoter_{idx}")
                    self.original_promoters[promoter_id] = {
                        'chr': row['chr'],
                        'start': row['start'],
                        'end': row['end'],
                        'strand': row.get('strand', '+'),
                        'name': promoter_id
                    }
                    
            logger.info(f"Loaded {len(self.original_promoters)} original promoter regions")
            
        except Exception as e:
            logger.error(f"Error reading original promoter file from {filename}: {e}")
            
    def find_nearest_tss(self, chr_name, start, end):
        """
        Find the nearest TSS to a given genomic region.
        
        Args:
            chr_name (str): Chromosome name
            start (int): Region start coordinate
            end (int): Region end coordinate
            
        Returns:
            dict: Nearest TSS information or None
        """
        nearest_tss = None
        min_distance = float('inf')
        
        region_center = (start + end) // 2
        
        for tss_region, tss_list in self.tss_sites.items():
            for tss in tss_list:
                if tss['chr'] == chr_name:
                    tss_center = (tss['start'] + tss['end']) // 2
                    distance = abs(region_center - tss_center)
                    
                    if distance < min_distance:
                        min_distance = distance
                        nearest_tss = tss.copy()
                        nearest_tss['distance'] = distance
                        
        return nearest_tss
        
    def check_histone_overlap(self, chr_name, start, end, mark_type="H3K4me3", min_overlap=0.5):
        """
        Check if a region overlaps with histone marks.
        
        Args:
            chr_name (str): Chromosome name
            start (int): Region start coordinate
            end (int): Region end coordinate
            mark_type (str): Type of histone mark to check
            min_overlap (float): Minimum fraction of overlap required
            
        Returns:
            bool: True if region has sufficient histone mark overlap
        """
        if mark_type not in self.histone_marks:
            return False
            
        region_length = end - start
        overlap_length = 0
        
        for mark in self.histone_marks[mark_type]:
            if mark['chr'] == chr_name:
                # Calculate overlap
                overlap_start = max(start, mark['start'])
                overlap_end = min(end, mark['end'])
                
                if overlap_start < overlap_end:
                    overlap_length += overlap_end - overlap_start
                    
        overlap_fraction = overlap_length / region_length if region_length > 0 else 0
        return overlap_fraction >= min_overlap
        
    def correct_promoter_coordinates(self, output_filename, min_histone_overlap=0.3):
        """
        Correct promoter coordinates based on TSS and histone marks.
        
        Args:
            output_filename (str): Path to output corrected promoter BED file
            min_histone_overlap (float): Minimum histone mark overlap required
        """
        corrected_promoters = []
        stats = {
            'total': len(self.original_promoters),
            'tss_corrected': 0,
            'histone_validated': 0,
            'no_tss_found': 0,
            'no_correction_needed': 0
        }
        
        logger.info("Starting promoter coordinate correction...")
        
        for promoter_id, promoter in self.original_promoters.items():
            chr_name = promoter['chr']
            start = promoter['start']
            end = promoter['end']
            strand = promoter['strand']
            
            # Find nearest TSS
            nearest_tss = self.find_nearest_tss(chr_name, start, end)
            
            if nearest_tss is None:
                logger.warning(f"No TSS found for promoter {promoter_id}")
                stats['no_tss_found'] += 1
                # Keep original coordinates if no TSS found
                corrected_promoters.append(promoter)
                continue
                
            # Calculate new promoter coordinates based on TSS
            tss_pos = (nearest_tss['start'] + nearest_tss['end']) // 2
            
            if strand == '-':
                # For negative strand, promoter is downstream of TSS
                new_start = tss_pos - self.downstream_bp
                new_end = tss_pos + self.upstream_bp
            else:
                # For positive strand, promoter is upstream of TSS
                new_start = tss_pos - self.upstream_bp
                new_end = tss_pos + self.downstream_bp
                
            # Ensure coordinates are positive
            new_start = max(0, new_start)
            
            # Check if correction is needed (significant change in coordinates)
            if abs(new_start - start) < 100 and abs(new_end - end) < 100:
                stats['no_correction_needed'] += 1
                corrected_promoters.append(promoter)
                continue
                
            # Validate with histone marks
            has_histone_marks = self.check_histone_overlap(
                chr_name, new_start, new_end, "H3K4me3", min_histone_overlap
            )
            
            if has_histone_marks:
                stats['histone_validated'] += 1
            else:
                # Try with alternative histone mark if available
                has_alt_marks = self.check_histone_overlap(
                    chr_name, new_start, new_end, "H3K27ac", min_histone_overlap
                )
                if has_alt_marks:
                    stats['histone_validated'] += 1
                    
            # Create corrected promoter entry
            corrected_promoter = {
                'chr': chr_name,
                'start': new_start,
                'end': new_end,
                'name': f"{promoter_id}_corrected",
                'score': nearest_tss.get('score', 0),
                'strand': strand,
                'original_start': start,
                'original_end': end,
                'tss_distance': nearest_tss['distance'],
                'histone_validated': has_histone_marks or has_alt_marks
            }
            
            corrected_promoters.append(corrected_promoter)
            stats['tss_corrected'] += 1
            
        # Write corrected promoters to output file
        self._write_corrected_bed(corrected_promoters, output_filename)
        
        # Print statistics
        logger.info("Promoter correction completed!")
        logger.info(f"Total promoters: {stats['total']}")
        logger.info(f"TSS-corrected: {stats['tss_corrected']}")
        logger.info(f"Histone-validated: {stats['histone_validated']}")
        logger.info(f"No TSS found: {stats['no_tss_found']}")
        logger.info(f"No correction needed: {stats['no_correction_needed']}")
        
        return corrected_promoters, stats
        
    def _write_corrected_bed(self, corrected_promoters, output_filename):
        """
        Write corrected promoter coordinates to BED file.
        
        Args:
            corrected_promoters (list): List of corrected promoter dictionaries
            output_filename (str): Path to output BED file
        """
        try:
            with open(output_filename, 'w') as f:
                # Write header comment
                f.write("# Corrected promoter coordinates\n")
                f.write("# Columns: chr, start, end, name, score, strand\n")
                
                for promoter in corrected_promoters:
                    line = f"{promoter['chr']}\t{promoter['start']}\t{promoter['end']}\t"
                    line += f"{promoter['name']}\t{promoter.get('score', 0)}\t{promoter['strand']}\n"
                    f.write(line)
                    
            logger.info(f"Corrected promoter coordinates written to {output_filename}")
            
        except Exception as e:
            logger.error(f"Error writing output file {output_filename}: {e}")


def main():
    """
    Main function to run the promoter correction tool.
    """
    parser = argparse.ArgumentParser(
        description="Correct promoter BED file coordinates using TSS and histone mark data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python promoter_corrector.py -p promoters.bed -g genes.gtf -t tss.bed --histone H3K4me3.bed -o corrected_promoters.bed
  
  python promoter_corrector.py -p promoters.bed -g genes.bed --gene-source refseq \\
                               -t cage_tss.bed --histone H3K4me3.bed --histone-type H3K4me3 \\
                               --upstream 2000 --downstream 500 -o output.bed
        """
    )
    
    # Required arguments
    parser.add_argument('-p', '--promoters', required=True,
                       help='Original promoter BED file to correct')
    parser.add_argument('-g', '--genes', required=True,
                       help='Gene coordinate file (GTF or BED format)')
    parser.add_argument('-t', '--tss', required=True,
                       help='CAGE TSS data file (BED format)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output corrected promoter BED file')
    
    # Optional arguments
    parser.add_argument('--histone', 
                       help='Histone mark data file (BED format)')
    parser.add_argument('--gene-source', choices=['gencode', 'refseq'], default='gencode',
                       help='Source of gene coordinates (default: gencode)')
    parser.add_argument('--histone-type', default='H3K4me3',
                       help='Type of histone mark (default: H3K4me3)')
    parser.add_argument('--upstream', type=int, default=1000,
                       help='Base pairs upstream of TSS (default: 1000)')
    parser.add_argument('--downstream', type=int, default=500,
                       help='Base pairs downstream of TSS (default: 500)')
    parser.add_argument('--min-overlap', type=float, default=0.3,
                       help='Minimum histone mark overlap fraction (default: 0.3)')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        
    # Initialize promoter corrector
    corrector = PromoterCorrector(args.upstream, args.downstream)
    
    # Load data files
    logger.info("Loading input data files...")
    
    # Load original promoters
    corrector.read_original_promoters(args.promoters)
    
    # Load gene coordinates
    corrector.read_gene_coordinates(args.genes, args.gene_source)
    
    # Load TSS data
    corrector.read_cage_tss_data(args.tss)
    
    # Load histone marks if provided
    if args.histone:
        corrector.read_histone_marks(args.histone, args.histone_type)
    
    # Perform correction
    logger.info("Performing promoter coordinate correction...")
    corrected_promoters, stats = corrector.correct_promoter_coordinates(
        args.output, args.min_overlap
    )
    
    logger.info("Promoter correction completed successfully!")


if __name__ == "__main__":
    main()