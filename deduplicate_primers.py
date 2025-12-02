#!/usr/bin/env python3
"""
Region-based Primer Deduplication
Removes primers that cover overlapping/nested regions, keeping only the best per region
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


class RegionDeduplicator:
    """Remove redundant primers covering the same regions"""
    
    def __init__(self, database_dir, primer_file, output_dir,
                 min_overlap_fraction=0.7,
                 min_primers_per_serotype=5,
                 max_mismatches=2):
        self.database_dir = Path(database_dir)
        self.primer_file = Path(primer_file)
        self.output_dir = Path(output_dir)
        self.min_overlap_fraction = min_overlap_fraction
        self.min_primers_per_serotype = min_primers_per_serotype
        self.max_mismatches = max_mismatches
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.reference_loci = {}
        self.primers = None
        self.primer_positions = {}  # primer -> serotype -> positions
        self.deduplicated_primers = []
        
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
    
    def map_primer_positions(self):
        """Map where each primer binds on each serotype"""
        logger.info("\nMapping primer positions...")
        
        for idx, primer_row in self.primers.iterrows():
            if (idx + 1) % 100 == 0:
                logger.info(f"  Mapped {idx + 1}/{len(self.primers)} primers...")
            
            primer_name = primer_row['name']
            left_seq = primer_row['left_seq']
            right_seq = primer_row['right_seq']
            right_seq_rc = str(Seq(right_seq).reverse_complement())
            
            self.primer_positions[primer_name] = {}
            
            for serotype, sequence in self.reference_loci.items():
                left_match = self._find_primer_in_sequence(left_seq, sequence, self.max_mismatches)
                right_match = self._find_primer_in_sequence(right_seq_rc, sequence, self.max_mismatches)
                
                if left_match and right_match and right_match['position'] > left_match['position']:
                    amplicon_start = left_match['position']
                    amplicon_end = right_match['position'] + len(right_seq)
                    amplicon_size = amplicon_end - amplicon_start
                    
                    self.primer_positions[primer_name][serotype] = {
                        'start': amplicon_start,
                        'end': amplicon_end,
                        'size': amplicon_size,
                        'mismatches': left_match['mismatches'] + right_match['mismatches']
                    }
        
        logger.info(f"Mapped positions for {len(self.primer_positions)} primers")
        return self.primer_positions
    
    def _intervals_overlap(self, interval1, interval2, min_fraction=0.7):
        """
        Check if two intervals overlap by at least min_fraction
        
        Returns True if overlap / shorter_interval >= min_fraction
        """
        start1, end1 = interval1
        start2, end2 = interval2
        
        # Calculate overlap
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        
        if overlap_end <= overlap_start:
            return False  # No overlap
        
        overlap_length = overlap_end - overlap_start
        
        # Calculate as fraction of shorter interval
        len1 = end1 - start1
        len2 = end2 - start2
        shorter_length = min(len1, len2)
        
        overlap_fraction = overlap_length / shorter_length if shorter_length > 0 else 0
        
        return overlap_fraction >= min_fraction
    
    def deduplicate_per_serotype(self):
        """
        For each serotype, remove primers that cover overlapping regions
        Keep the "best" primer per region (based on discriminatory score or mismatches)
        """
        logger.info("\n" + "=" * 80)
        logger.info("DEDUPLICATING OVERLAPPING PRIMERS")
        logger.info("=" * 80)
        logger.info(f"Overlap threshold: {self.min_overlap_fraction * 100:.0f}% of shorter amplicon")
        
        # Track which primers to keep
        primers_to_keep = set()
        
        # Process each serotype
        for serotype in sorted(self.reference_loci.keys()):
            logger.info(f"\nProcessing {serotype}...")
            
            # Get all primers that amplify this serotype
            primers_for_serotype = []
            for primer_name, positions in self.primer_positions.items():
                if serotype in positions:
                    pos = positions[serotype]
                    primers_for_serotype.append({
                        'name': primer_name,
                        'start': pos['start'],
                        'end': pos['end'],
                        'size': pos['size'],
                        'mismatches': pos['mismatches']
                    })
            
            if not primers_for_serotype:
                logger.info(f"  No primers found for {serotype}")
                continue
            
            logger.info(f"  Found {len(primers_for_serotype)} primers")
            
            # Sort by position
            primers_for_serotype.sort(key=lambda x: x['start'])
            
            # Greedy deduplication: keep best non-overlapping primers
            selected_primers = []
            
            for primer in primers_for_serotype:
                # Check if this primer overlaps with any already selected
                overlaps = False
                
                for selected in selected_primers:
                    if self._intervals_overlap(
                        (primer['start'], primer['end']),
                        (selected['start'], selected['end']),
                        self.min_overlap_fraction
                    ):
                        overlaps = True
                        break
                
                if not overlaps:
                    selected_primers.append(primer)
                    primers_to_keep.add(primer['name'])
            
            logger.info(f"  After deduplication: {len(selected_primers)} primers (removed {len(primers_for_serotype) - len(selected_primers)})")
            
            # Check if we have minimum coverage
            if len(selected_primers) < self.min_primers_per_serotype:
                logger.warning(f"  ⚠ {serotype} has only {len(selected_primers)} non-overlapping primers (min: {self.min_primers_per_serotype})")
                
                # Add back some overlapping primers if needed
                remaining = [p for p in primers_for_serotype if p['name'] not in primers_to_keep]
                remaining.sort(key=lambda x: x['mismatches'])  # Best matches first
                
                for primer in remaining:
                    if len(selected_primers) >= self.min_primers_per_serotype:
                        break
                    selected_primers.append(primer)
                    primers_to_keep.add(primer['name'])
                
                logger.info(f"  Added {len(selected_primers) - len([p for p in primers_for_serotype if p['name'] in primers_to_keep and p in selected_primers])} overlapping primers to reach minimum")
        
        self.deduplicated_primers = list(primers_to_keep)
        
        logger.info("\n" + "=" * 80)
        logger.info("DEDUPLICATION SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Original primers: {len(self.primers)}")
        logger.info(f"Deduplicated primers: {len(self.deduplicated_primers)}")
        logger.info(f"Removed: {len(self.primers) - len(self.deduplicated_primers)} ({(1 - len(self.deduplicated_primers)/len(self.primers))*100:.1f}%)")
        
        return self.deduplicated_primers
    
    def export_deduplicated_primers(self):
        """Export deduplicated primers"""
        logger.info("\nExporting deduplicated primers...")
        
        dedup_df = self.primers[self.primers['name'].isin(self.deduplicated_primers)]
        
        output_file = self.output_dir / "deduplicated_primers.csv"
        dedup_df.to_csv(output_file, index=False)
        logger.info(f"  Saved: {output_file}")
        
        # Export by pool
        for pool_name in ['pool1', 'pool2']:
            pool_df = dedup_df[dedup_df['pool'] == pool_name]
            if len(pool_df) > 0:
                pool_file = self.output_dir / f"{pool_name}_deduplicated.csv"
                pool_df.to_csv(pool_file, index=False)
                logger.info(f"  {pool_name}: {pool_file} ({len(pool_df)} primers)")
        
        return output_file
    
    def generate_coverage_report(self):
        """Generate report showing coverage per serotype"""
        report_file = self.output_dir / "deduplication_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("REGION DEDUPLICATION REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("PARAMETERS:\n")
            f.write(f"  Overlap threshold: {self.min_overlap_fraction * 100:.0f}%\n")
            f.write(f"  Min primers per serotype: {self.min_primers_per_serotype}\n\n")
            
            f.write("RESULTS:\n")
            f.write(f"  Input primers: {len(self.primers)}\n")
            f.write(f"  Output primers: {len(self.deduplicated_primers)}\n")
            f.write(f"  Reduction: {(1 - len(self.deduplicated_primers)/len(self.primers))*100:.1f}%\n\n")
            
            f.write("COVERAGE PER SEROTYPE:\n")
            
            for serotype in sorted(self.reference_loci.keys()):
                # Count primers for this serotype
                count = sum(1 for name in self.deduplicated_primers 
                           if serotype in self.primer_positions.get(name, {}))
                
                status = "OK" if count >= self.min_primers_per_serotype else "LOW"
                f.write(f"  {serotype:30s}: {count:3d} primers [{status}]\n")
            
            f.write("\n" + "=" * 80 + "\n")
        
        logger.info(f"Report saved: {report_file}")
        return report_file


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Remove redundant primers covering overlapping regions'
    )
    
    parser.add_argument('--database', required=True)
    parser.add_argument('--primers', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--overlap-threshold', type=float, default=0.7,
                       help='Min overlap fraction to consider redundant (default: 0.7)')
    parser.add_argument('--min-per-serotype', type=int, default=5,
                       help='Min primers per serotype (default: 5)')
    parser.add_argument('--max-mismatches', type=int, default=2,
                       help='Max primer mismatches (default: 2)')
    parser.add_argument('--debug', action='store_true')
    
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    deduplicator = RegionDeduplicator(
        database_dir=args.database,
        primer_file=args.primers,
        output_dir=args.output,
        min_overlap_fraction=args.overlap_threshold,
        min_primers_per_serotype=args.min_per_serotype,
        max_mismatches=args.max_mismatches
    )
    
    try:
        logger.info("=" * 80)
        logger.info("REGION-BASED PRIMER DEDUPLICATION")
        logger.info("=" * 80)
        
        deduplicator.load_reference_loci()
        deduplicator.load_primers()
        deduplicator.map_primer_positions()
        deduplicator.deduplicate_per_serotype()
        deduplicator.export_deduplicated_primers()
        deduplicator.generate_coverage_report()
        
        logger.info("=" * 80)
        logger.info("COMPLETE!")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Failed: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
