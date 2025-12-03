#!/usr/bin/env python3
"""
Deduplication using Actual Primer Mapping
Maps primers using in silico PCR, then deduplicates based on actual amplicon positions
"""

import os
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class TrueDeduplicator:
    """Deduplicate using actual primer mapping, not metadata"""
    
    def __init__(self, database_dir, primer_file, output_dir,
                 min_gap=100, min_primers_per_serotype=5, max_mismatches=2):
        self.database_dir = Path(database_dir)
        self.primer_file = Path(primer_file)
        self.output_dir = Path(output_dir)
        self.min_gap = min_gap
        self.min_primers_per_serotype = min_primers_per_serotype
        self.max_mismatches = max_mismatches
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.reference_loci = {}
        self.primers = None
        self.actual_amplicons = {}  # primer -> serotype -> {start, end, ...}
        self.selected_primers = set()
        
    def load_reference_loci(self):
        """Load reference loci"""
        logger.info("Loading reference loci...")
        
        root_ref = self.database_dir / "reference.fasta"
        if root_ref.exists():
            for record in SeqIO.parse(root_ref, "fasta"):
                self.reference_loci[record.id] = str(record.seq).upper()
        
        for subdir in self.database_dir.iterdir():
            if subdir.is_dir():
                ref_file = subdir / "reference.fasta"
                if ref_file.exists():
                    for record in SeqIO.parse(ref_file, "fasta"):
                        serotype = subdir.name
                        if serotype not in self.reference_loci:
                            self.reference_loci[serotype] = str(record.seq).upper()
        
        logger.info(f"Loaded {len(self.reference_loci)} reference loci")
        return self.reference_loci
    
    def load_primers(self):
        """Load primers"""
        logger.info(f"Loading primers from {self.primer_file}")
        
        if self.primer_file.suffix == '.tsv':
            self.primers = pd.read_csv(self.primer_file, sep='\t')
        else:
            self.primers = pd.read_csv(self.primer_file)
        
        logger.info(f"Loaded {len(self.primers)} primers")
        return self.primers
    
    def _find_primer_in_sequence(self, primer_seq, sequence, max_mismatches):
        """Find primer binding site"""
        primer_len = len(primer_seq)
        best_match = None
        
        for i in range(len(sequence) - primer_len + 1):
            target_substr = sequence[i:i + primer_len]
            mismatches = sum(1 for a, b in zip(primer_seq, target_substr) if a != b)
            
            if mismatches <= max_mismatches:
                if best_match is None or mismatches < best_match['mismatches']:
                    best_match = {'position': i, 'mismatches': mismatches}
        
        return best_match
    
    def map_all_primers(self):
        """
        Perform in silico PCR for ALL primers on ALL serotypes
        This is the ground truth of what actually amplifies
        """
        logger.info("\nPerforming in silico PCR for all primers on all serotypes...")
        logger.info("This may take a while...")
        
        for idx, primer_row in self.primers.iterrows():
            if (idx + 1) % 50 == 0:
                logger.info(f"  Mapped {idx + 1}/{len(self.primers)} primers...")
            
            primer_name = primer_row['name']
            left_seq = primer_row['left_seq']
            right_seq = primer_row['right_seq']
            right_seq_rc = str(Seq(right_seq).reverse_complement())
            
            self.actual_amplicons[primer_name] = {}
            
            for serotype, sequence in self.reference_loci.items():
                left_match = self._find_primer_in_sequence(left_seq, sequence, self.max_mismatches)
                right_match = self._find_primer_in_sequence(right_seq_rc, sequence, self.max_mismatches)
                
                if left_match and right_match and right_match['position'] > left_match['position']:
                    amplicon_start = left_match['position']
                    amplicon_end = right_match['position'] + len(right_seq)
                    amplicon_size = amplicon_end - amplicon_start
                    
                    # Only keep reasonable amplicon sizes
                    if 200 <= amplicon_size <= 1000:
                        self.actual_amplicons[primer_name][serotype] = {
                            'start': amplicon_start,
                            'end': amplicon_end,
                            'size': amplicon_size,
                            'mismatches': left_match['mismatches'] + right_match['mismatches']
                        }
        
        # Count how many primers amplify each serotype
        logger.info("\nPrimers per serotype (before deduplication):")
        for serotype in sorted(self.reference_loci.keys())[:20]:  # Show first 20
            count = sum(1 for p in self.actual_amplicons.values() if serotype in p)
            logger.info(f"  {serotype:20s}: {count:3d} primers")
        logger.info("  ...")
        
        return self.actual_amplicons
    
    def _intervals_overlap(self, interval1, interval2, min_gap):
        """Check if two intervals overlap"""
        start1, end1 = interval1
        start2, end2 = interval2
        
        if start2 <= end1 and start1 <= end2:
            return True
        
        if start2 > end1:
            gap = start2 - end1
        else:
            gap = start1 - end2
        
        return gap < min_gap
    
    def deduplicate_globally(self):
        """
        Global deduplication across all serotypes
        Select primers that create non-overlapping amplicons in EVERY serotype they amplify
        """
        logger.info("\n" + "=" * 80)
        logger.info("GLOBAL DEDUPLICATION")
        logger.info("=" * 80)
        logger.info(f"Strategy: Greedy selection ensuring no overlaps in ANY serotype")
        
        # Build per-serotype primer lists
        serotype_primers = defaultdict(list)
        for primer_name, serotypes_dict in self.actual_amplicons.items():
            for serotype, pos in serotypes_dict.items():
                serotype_primers[serotype].append({
                    'name': primer_name,
                    'start': pos['start'],
                    'end': pos['end'],
                    'mismatches': pos['mismatches']
                })
        
        # Greedy selection
        selected_primers = set()
        serotype_coverage = defaultdict(int)
        serotype_covered_regions = defaultdict(list)  # serotype -> list of (start, end)
        
        # Sort all primers by how many serotypes they cover (prefer universal)
        primer_serotype_counts = [(p, len(amps)) for p, amps in self.actual_amplicons.items()]
        primer_serotype_counts.sort(key=lambda x: x[1], reverse=True)
        
        logger.info("\nPhase 1: Greedy selection...")
        
        for primer_name, num_serotypes in primer_serotype_counts:
            # Check if this primer overlaps with any already-selected primer in ANY serotype
            causes_overlap = False
            
            for serotype, pos in self.actual_amplicons[primer_name].items():
                # Check against already covered regions in this serotype
                for covered_start, covered_end in serotype_covered_regions[serotype]:
                    if self._intervals_overlap(
                        (pos['start'], pos['end']),
                        (covered_start, covered_end),
                        self.min_gap
                    ):
                        causes_overlap = True
                        break
                
                if causes_overlap:
                    break
            
            if not causes_overlap:
                # Add this primer
                selected_primers.add(primer_name)
                
                # Mark regions as covered
                for serotype, pos in self.actual_amplicons[primer_name].items():
                    serotype_covered_regions[serotype].append((pos['start'], pos['end']))
                    serotype_coverage[serotype] += 1
        
        logger.info(f"Phase 1 complete: {len(selected_primers)} primers selected")
        
        # Check coverage
        logger.info("\nCoverage per serotype:")
        low_coverage_serotypes = []
        for serotype in sorted(self.reference_loci.keys()):
            count = serotype_coverage.get(serotype, 0)
            status = "✓" if count >= self.min_primers_per_serotype else "⚠"
            if count < self.min_primers_per_serotype:
                low_coverage_serotypes.append(serotype)
            
            if count > 0 or serotype in list(self.reference_loci.keys())[:20]:
                logger.info(f"  {status} {serotype:20s}: {count:3d} primers")
        
        if low_coverage_serotypes:
            logger.warning(f"\n{len(low_coverage_serotypes)} serotypes have insufficient coverage")
            logger.warning("Consider lowering --min-gap to allow closer amplicons")
        
        self.selected_primers = selected_primers
        
        logger.info("\n" + "=" * 80)
        logger.info("DEDUPLICATION SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Original primers: {len(self.primers)}")
        logger.info(f"Selected primers: {len(self.selected_primers)}")
        logger.info(f"Removed: {len(self.primers) - len(self.selected_primers)} "
                   f"({(1 - len(self.selected_primers)/len(self.primers))*100:.1f}%)")
        
        return list(self.selected_primers)
    
    def export_primers(self):
        """Export deduplicated primers"""
        logger.info("\nExporting deduplicated primers...")
        
        dedup_df = self.primers[self.primers['name'].isin(self.selected_primers)]
        
        output_file = self.output_dir / "truly_deduplicated_primers.csv"
        dedup_df.to_csv(output_file, index=False)
        logger.info(f"  Saved: {output_file}")
        
        # Export by pool
        for pool_name in ['pool1', 'pool2']:
            if 'pool' in dedup_df.columns:
                pool_df = dedup_df[dedup_df['pool'] == pool_name]
                if len(pool_df) > 0:
                    pool_file = self.output_dir / f"{pool_name}_truly_deduplicated.csv"
                    pool_df.to_csv(pool_file, index=False)
                    logger.info(f"  {pool_name}: {pool_file} ({len(pool_df)} primers)")
        
        return output_file
    
    def generate_report(self):
        """Generate report"""
        report_file = self.output_dir / "true_deduplication_report.txt"
        
        # Calculate coverage per serotype
        serotype_counts = defaultdict(int)
        for primer_name in self.selected_primers:
            for serotype in self.actual_amplicons[primer_name].keys():
                serotype_counts[serotype] += 1
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("TRUE DEDUPLICATION REPORT\n")
            f.write("(Using actual in silico PCR mapping)\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("METHOD:\n")
            f.write("  1. Mapped ALL primers to ALL serotypes using in silico PCR\n")
            f.write("  2. Greedy selection ensuring no overlaps in ANY serotype\n")
            f.write("  3. Result: Globally non-overlapping primer set\n\n")
            
            f.write("PARAMETERS:\n")
            f.write(f"  Minimum gap: {self.min_gap}bp\n")
            f.write(f"  Max mismatches: {self.max_mismatches}\n")
            f.write(f"  Min primers per serotype: {self.min_primers_per_serotype}\n\n")
            
            f.write("RESULTS:\n")
            f.write(f"  Input primers: {len(self.primers)}\n")
            f.write(f"  Output primers: {len(self.selected_primers)}\n")
            f.write(f"  Reduction: {(1 - len(self.selected_primers)/len(self.primers))*100:.1f}%\n\n")
            
            f.write("NON-OVERLAPPING PRIMERS PER SEROTYPE:\n")
            for serotype in sorted(self.reference_loci.keys()):
                count = serotype_counts.get(serotype, 0)
                status = "OK" if count >= self.min_primers_per_serotype else "LOW"
                f.write(f"  {serotype:30s}: {count:3d} primers [{status}]\n")
            
            f.write("\n" + "=" * 80 + "\n")
        
        logger.info(f"Report saved: {report_file}")
        return report_file


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='True deduplication using actual primer mapping'
    )
    
    parser.add_argument('--database', required=True)
    parser.add_argument('--primers', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--min-gap', type=int, default=100)
    parser.add_argument('--min-per-serotype', type=int, default=5)
    parser.add_argument('--max-mismatches', type=int, default=2)
    parser.add_argument('--debug', action='store_true')
    
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    deduplicator = TrueDeduplicator(
        database_dir=args.database,
        primer_file=args.primers,
        output_dir=args.output,
        min_gap=args.min_gap,
        min_primers_per_serotype=args.min_per_serotype,
        max_mismatches=args.max_mismatches
    )
    
    try:
        logger.info("=" * 80)
        logger.info("TRUE DEDUPLICATION (USING ACTUAL PRIMER MAPPING)")
        logger.info("=" * 80)
        
        deduplicator.load_reference_loci()
        deduplicator.load_primers()
        deduplicator.map_all_primers()  # Key difference: actual mapping!
        deduplicator.deduplicate_globally()
        deduplicator.export_primers()
        deduplicator.generate_report()
        
        logger.info("=" * 80)
        logger.info("COMPLETE!")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Failed: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
