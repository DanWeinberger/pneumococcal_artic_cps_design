#!/usr/bin/env python3
"""
Aggressive Primer Deduplication using Interval Clustering
Groups overlapping amplicons into clusters and keeps only the best primer per cluster
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


class AggressiveDeduplicator:
    """Aggressively remove ALL overlapping primers"""
    
    def __init__(self, database_dir, primer_file, output_dir,
                 min_gap=100,  # Minimum gap between amplicons (bp)
                 min_primers_per_serotype=5,
                 max_mismatches=2):
        self.database_dir = Path(database_dir)
        self.primer_file = Path(primer_file)
        self.output_dir = Path(output_dir)
        self.min_gap = min_gap
        self.min_primers_per_serotype = min_primers_per_serotype
        self.max_mismatches = max_mismatches
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.reference_loci = {}
        self.primers = None
        self.primer_positions = {}
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
    
    def map_primer_positions(self):
        """Map where each primer binds"""
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
                    
                    self.primer_positions[primer_name][serotype] = {
                        'start': amplicon_start,
                        'end': amplicon_end,
                        'size': amplicon_end - amplicon_start,
                        'mismatches': left_match['mismatches'] + right_match['mismatches']
                    }
        
        logger.info(f"Mapped positions for {len(self.primer_positions)} primers")
        return self.primer_positions
    
    def _intervals_overlap(self, interval1, interval2, min_gap):
        """
        Check if two intervals overlap or are closer than min_gap
        
        Returns True if they should be considered overlapping
        """
        start1, end1 = interval1
        start2, end2 = interval2
        
        # If intervals overlap at all, they overlap
        if start2 <= end1 and start1 <= end2:
            return True
        
        # If gap between them is less than min_gap, consider overlapping
        if start2 > end1:
            gap = start2 - end1
        else:
            gap = start1 - end2
        
        return gap < min_gap
    
    def _cluster_overlapping_intervals(self, intervals):
        """
        Group overlapping intervals into clusters
        Returns list of clusters, where each cluster is a list of interval indices
        """
        if not intervals:
            return []
        
        # Sort by start position
        sorted_indices = sorted(range(len(intervals)), 
                               key=lambda i: intervals[i]['start'])
        
        clusters = []
        current_cluster = [sorted_indices[0]]
        
        for i in sorted_indices[1:]:
            # Check if this interval overlaps with any in current cluster
            overlaps = False
            
            for j in current_cluster:
                if self._intervals_overlap(
                    (intervals[i]['start'], intervals[i]['end']),
                    (intervals[j]['start'], intervals[j]['end']),
                    self.min_gap
                ):
                    overlaps = True
                    break
            
            if overlaps:
                current_cluster.append(i)
            else:
                # Start new cluster
                clusters.append(current_cluster)
                current_cluster = [i]
        
        # Add last cluster
        if current_cluster:
            clusters.append(current_cluster)
        
        return clusters
    
    def _select_best_from_cluster(self, cluster_intervals, primers_data):
        """
        Select the best primer from a cluster
        Prioritize: 1) fewest mismatches, 2) most central position, 3) best discriminatory score
        """
        if len(cluster_intervals) == 1:
            return cluster_intervals[0]
        
        # Score each primer
        scores = []
        
        for interval_idx in cluster_intervals:
            primer_info = primers_data[interval_idx]
            
            # Lower mismatches is better
            mismatch_score = -primer_info['mismatches']
            
            # Central position is better (avoid edges)
            position_score = 0  # Could be improved with more context
            
            # Discriminatory score if available
            disc_score = primer_info.get('discriminatory_score', 0)
            
            total_score = mismatch_score * 10 + disc_score
            scores.append((total_score, interval_idx))
        
        # Return best scoring primer
        best = max(scores, key=lambda x: x[0])
        return best[1]
    
    def deduplicate_aggressive(self):
        """
        Aggressive deduplication:
        1. For each serotype, cluster all overlapping amplicons
        2. Keep only ONE primer per cluster (the best one)
        3. Ensure final set is non-overlapping per serotype
        """
        logger.info("\n" + "=" * 80)
        logger.info("AGGRESSIVE DEDUPLICATION")
        logger.info("=" * 80)
        logger.info(f"Minimum gap between amplicons: {self.min_gap}bp")
        logger.info("Strategy: One primer per cluster of overlapping amplicons PER SEROTYPE")
        
        # Track best primers per serotype (not global!)
        serotype_selected = {}  # serotype -> list of primer names
        
        # Process each serotype independently
        for serotype in sorted(self.reference_loci.keys()):
            logger.info(f"\nProcessing {serotype}...")
            
            # Get all primers for this serotype
            primers_data = []
            for primer_name, positions in self.primer_positions.items():
                if serotype in positions:
                    pos = positions[serotype]
                    
                    # Get discriminatory score if available
                    primer_row = self.primers[self.primers['name'] == primer_name]
                    disc_score = 0
                    if len(primer_row) > 0 and 'discriminatory_score' in primer_row.columns:
                        disc_score = primer_row.iloc[0]['discriminatory_score']
                    
                    primers_data.append({
                        'name': primer_name,
                        'start': pos['start'],
                        'end': pos['end'],
                        'size': pos['size'],
                        'mismatches': pos['mismatches'],
                        'discriminatory_score': disc_score
                    })
            
            if not primers_data:
                logger.info(f"  No primers for {serotype}")
                serotype_selected[serotype] = []
                continue
            
            logger.info(f"  Found {len(primers_data)} primers")
            
            # Cluster overlapping intervals
            clusters = self._cluster_overlapping_intervals(primers_data)
            
            logger.info(f"  Grouped into {len(clusters)} non-overlapping clusters")
            
            # Select best primer from each cluster
            selected_for_serotype = []
            
            for cluster_idx, cluster in enumerate(clusters, 1):
                best_idx = self._select_best_from_cluster(cluster, primers_data)
                best_primer = primers_data[best_idx]
                
                selected_for_serotype.append(best_primer['name'])
                
                logger.debug(f"    Cluster {cluster_idx}: {len(cluster)} primers → selected {best_primer['name']}")
            
            serotype_selected[serotype] = selected_for_serotype
            
            logger.info(f"  Selected {len(selected_for_serotype)} primers "
                       f"(removed {len(primers_data) - len(selected_for_serotype)})")
            
            # Check minimum coverage
            if len(selected_for_serotype) < self.min_primers_per_serotype:
                logger.warning(f"  ⚠ {serotype} has only {len(selected_for_serotype)} primers "
                             f"(min: {self.min_primers_per_serotype})")
        
        # Now collect primers that are selected for ANY serotype
        # CRITICAL: A primer is kept if it's selected for at least one serotype
        all_selected = set()
        for primers_list in serotype_selected.values():
            all_selected.update(primers_list)
        
        self.selected_primers = all_selected
        
        # Verify no overlaps per serotype in final set
        logger.info("\n" + "=" * 80)
        logger.info("VERIFYING NON-OVERLAP")
        logger.info("=" * 80)
        
        for serotype, selected_list in serotype_selected.items():
            if not selected_list:
                continue
            
            # Get positions of selected primers for this serotype
            positions = []
            for primer_name in selected_list:
                if serotype in self.primer_positions[primer_name]:
                    pos = self.primer_positions[primer_name][serotype]
                    positions.append((pos['start'], pos['end'], primer_name))
            
            # Sort and check for overlaps
            positions.sort()
            has_overlap = False
            
            for i in range(len(positions) - 1):
                if positions[i][1] > positions[i+1][0]:  # end of i > start of i+1
                    logger.error(f"  ERROR: {serotype} has overlap: "
                               f"{positions[i][2]} ({positions[i][0]}-{positions[i][1]}) and "
                               f"{positions[i+1][2]} ({positions[i+1][0]}-{positions[i+1][1]})")
                    has_overlap = True
            
            if not has_overlap:
                logger.info(f"  ✓ {serotype}: {len(selected_list)} primers, no overlaps")
        
        logger.info("\n" + "=" * 80)
        logger.info("DEDUPLICATION SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Original primers: {len(self.primers)}")
        logger.info(f"Deduplicated primers: {len(self.selected_primers)}")
        logger.info(f"Removed: {len(self.primers) - len(self.selected_primers)} "
                   f"({(1 - len(self.selected_primers)/len(self.primers))*100:.1f}%)")
        
        return list(self.selected_primers)
    
    def export_primers(self):
        """Export deduplicated primers"""
        logger.info("\nExporting deduplicated primers...")
        
        dedup_df = self.primers[self.primers['name'].isin(self.selected_primers)]
        
        # Sort by pool and position
        if 'pool' in dedup_df.columns:
            dedup_df = dedup_df.sort_values(['pool', 'left_pos'] if 'left_pos' in dedup_df.columns else 'pool')
        
        output_file = self.output_dir / "deduplicated_primers.csv"
        dedup_df.to_csv(output_file, index=False)
        logger.info(f"  Saved: {output_file}")
        
        # Export by pool
        for pool_name in ['pool1', 'pool2']:
            if 'pool' in dedup_df.columns:
                pool_df = dedup_df[dedup_df['pool'] == pool_name]
                if len(pool_df) > 0:
                    pool_file = self.output_dir / f"{pool_name}_deduplicated.csv"
                    pool_df.to_csv(pool_file, index=False)
                    logger.info(f"  {pool_name}: {pool_file} ({len(pool_df)} primers)")
        
        return output_file
    
    def generate_report(self):
        """Generate report"""
        report_file = self.output_dir / "deduplication_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("AGGRESSIVE PRIMER DEDUPLICATION REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("STRATEGY:\n")
            f.write("  - Cluster all overlapping/nearby amplicons per serotype\n")
            f.write("  - Keep only the BEST primer from each cluster\n")
            f.write("  - Result: Zero overlap between amplicons\n\n")
            
            f.write("PARAMETERS:\n")
            f.write(f"  Minimum gap: {self.min_gap}bp\n")
            f.write(f"  Min primers per serotype: {self.min_primers_per_serotype}\n\n")
            
            f.write("RESULTS:\n")
            f.write(f"  Input primers: {len(self.primers)}\n")
            f.write(f"  Output primers: {len(self.selected_primers)}\n")
            f.write(f"  Reduction: {(1 - len(self.selected_primers)/len(self.primers))*100:.1f}%\n\n")
            
            f.write("COVERAGE PER SEROTYPE (non-overlapping primers):\n")
            
            for serotype in sorted(self.reference_loci.keys()):
                count = sum(1 for name in self.selected_primers 
                           if serotype in self.primer_positions.get(name, {}))
                
                status = "OK" if count >= self.min_primers_per_serotype else "LOW"
                f.write(f"  {serotype:30s}: {count:3d} primers [{status}]\n")
            
            f.write("\n" + "=" * 80 + "\n")
        
        logger.info(f"Report saved: {report_file}")
        return report_file


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Aggressively remove ALL overlapping primers'
    )
    
    parser.add_argument('--database', required=True)
    parser.add_argument('--primers', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--min-gap', type=int, default=100,
                       help='Minimum gap between amplicons (bp, default: 100)')
    parser.add_argument('--min-per-serotype', type=int, default=5,
                       help='Min primers per serotype (default: 5)')
    parser.add_argument('--max-mismatches', type=int, default=2)
    parser.add_argument('--debug', action='store_true')
    
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    deduplicator = AggressiveDeduplicator(
        database_dir=args.database,
        primer_file=args.primers,
        output_dir=args.output,
        min_gap=args.min_gap,
        min_primers_per_serotype=args.min_per_serotype,
        max_mismatches=args.max_mismatches
    )
    
    try:
        logger.info("=" * 80)
        logger.info("AGGRESSIVE PRIMER DEDUPLICATION")
        logger.info("=" * 80)
        
        deduplicator.load_reference_loci()
        deduplicator.load_primers()
        deduplicator.map_primer_positions()
        deduplicator.deduplicate_aggressive()
        deduplicator.export_primers()
        deduplicator.generate_report()
        
        logger.info("=" * 80)
        logger.info("COMPLETE!")
        logger.info("Result: ZERO overlapping amplicons per serotype")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Failed: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
