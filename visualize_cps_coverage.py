#!/usr/bin/env python3
"""
CPS Locus Coverage Visualizer
Maps amplicons to reference CPS loci and visualizes physical coverage
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
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class CPSCoverageVisualizer:
    """Visualize amplicon coverage across CPS loci"""
    
    def __init__(self, database_dir, primer_file, output_dir, max_mismatches=2):
        self.database_dir = Path(database_dir)
        self.primer_file = Path(primer_file)
        self.output_dir = Path(output_dir)
        self.max_mismatches = max_mismatches
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "plots").mkdir(exist_ok=True)
        (self.output_dir / "coverage_reports").mkdir(exist_ok=True)
        
        self.primers = None
        self.reference_loci = {}
        self.amplicon_coverage = defaultdict(list)
        
    def load_primers(self):
        """Load primer sequences"""
        logger.info(f"Loading primers from {self.primer_file}")
        
        if self.primer_file.suffix == '.tsv':
            self.primers = pd.read_csv(self.primer_file, sep='\t')
        else:
            self.primers = pd.read_csv(self.primer_file)
        
        logger.info(f"Loaded {len(self.primers)} primer pairs")
        return self.primers
    
    def load_reference_loci(self):
        """Load reference CPS loci (whole cassettes)"""
        logger.info("Loading reference CPS loci...")
        
        # Load root reference.fasta (whole cassettes)
        root_ref = self.database_dir / "reference.fasta"
        if root_ref.exists():
            for record in SeqIO.parse(root_ref, "fasta"):
                serotype = record.id
                self.reference_loci[serotype] = {
                    'id': record.id,
                    'sequence': str(record.seq).upper(),
                    'length': len(record.seq),
                    'source': 'root'
                }
        
        # Load subdirectory reference.fasta files
        for subdir in self.database_dir.iterdir():
            if subdir.is_dir():
                ref_file = subdir / "reference.fasta"
                if ref_file.exists():
                    for record in SeqIO.parse(ref_file, "fasta"):
                        serotype = subdir.name
                        # Only add if not already present (root takes precedence)
                        if serotype not in self.reference_loci:
                            self.reference_loci[serotype] = {
                                'id': record.id,
                                'sequence': str(record.seq).upper(),
                                'length': len(record.seq),
                                'source': subdir.name
                            }
        
        logger.info(f"Loaded {len(self.reference_loci)} reference CPS loci")
        
        # Show locus sizes
        logger.info("\nReference locus sizes:")
        for serotype, data in sorted(self.reference_loci.items())[:10]:
            logger.info(f"  {serotype}: {data['length']:,} bp")
        
        return self.reference_loci
    
    def _find_primer_in_sequence(self, primer_seq, sequence, max_mismatches):
        """Find primer binding site in sequence"""
        primer_len = len(primer_seq)
        best_match = None
        
        for i in range(len(sequence) - primer_len + 1):
            target_substr = sequence[i:i + primer_len]
            mismatches = sum(1 for a, b in zip(primer_seq, target_substr) if a != b)
            
            if mismatches <= max_mismatches:
                if best_match is None or mismatches < best_match['mismatches']:
                    best_match = {
                        'position': i,
                        'mismatches': mismatches,
                        'sequence': target_substr
                    }
        
        return best_match
    
    def map_amplicons_to_loci(self):
        """Map each amplicon to reference CPS loci"""
        logger.info("Mapping amplicons to reference CPS loci...")
        
        coverage_data = []
        
        for serotype, locus_data in self.reference_loci.items():
            logger.info(f"Processing {serotype} ({locus_data['length']:,} bp)...")
            
            sequence = locus_data['sequence']
            amplicons = []
            
            # Try each primer pair on this locus
            for _, primer in self.primers.iterrows():
                left_seq = primer['left_seq']
                right_seq = primer['right_seq']
                right_seq_rc = str(Seq(right_seq).reverse_complement())
                
                # Find binding sites
                left_match = self._find_primer_in_sequence(left_seq, sequence, self.max_mismatches)
                right_match = self._find_primer_in_sequence(right_seq_rc, sequence, self.max_mismatches)
                
                # Check if both primers bind and create valid amplicon
                if left_match and right_match:
                    if right_match['position'] > left_match['position']:
                        amplicon_start = left_match['position']
                        amplicon_end = right_match['position'] + len(right_seq)
                        amplicon_size = amplicon_end - amplicon_start
                        
                        # Check if size is reasonable
                        expected_size = primer['product_size']
                        if 0.3 * expected_size <= amplicon_size <= 3.0 * expected_size:
                            amplicons.append({
                                'primer_name': primer['name'],
                                'gene': primer['gene'],
                                'start': amplicon_start,
                                'end': amplicon_end,
                                'size': amplicon_size,
                                'left_mismatches': left_match['mismatches'],
                                'right_mismatches': right_match['mismatches'],
                                'total_mismatches': left_match['mismatches'] + right_match['mismatches']
                            })
            
            self.amplicon_coverage[serotype] = amplicons
            
            # Calculate coverage statistics
            if amplicons:
                covered_positions = set()
                for amp in amplicons:
                    covered_positions.update(range(amp['start'], amp['end']))
                
                coverage_percent = (len(covered_positions) / locus_data['length']) * 100
                
                coverage_data.append({
                    'serotype': serotype,
                    'locus_length': locus_data['length'],
                    'num_amplicons': len(amplicons),
                    'covered_bases': len(covered_positions),
                    'coverage_percent': coverage_percent
                })
                
                logger.info(f"  {len(amplicons)} amplicons, {coverage_percent:.1f}% coverage")
            else:
                coverage_data.append({
                    'serotype': serotype,
                    'locus_length': locus_data['length'],
                    'num_amplicons': 0,
                    'covered_bases': 0,
                    'coverage_percent': 0.0
                })
                logger.warning(f"  No amplicons found!")
        
        self.coverage_summary = pd.DataFrame(coverage_data)
        
        return self.coverage_summary
    
    def plot_locus_coverage(self, serotype, figsize=(16, 4)):
        """Plot amplicon coverage for a single serotype"""
        
        if serotype not in self.reference_loci:
            logger.error(f"Serotype {serotype} not found in reference loci")
            return None
        
        locus_data = self.reference_loci[serotype]
        amplicons = self.amplicon_coverage[serotype]
        locus_length = locus_data['length']
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Draw reference locus as gray bar
        ax.add_patch(Rectangle((0, 0), locus_length, 0.5, 
                               facecolor='lightgray', edgecolor='black', linewidth=1))
        
        # Draw amplicons
        colors = plt.cm.tab20(np.linspace(0, 1, 20))
        
        for i, amp in enumerate(amplicons):
            color = colors[i % 20]
            y_offset = 0.6 + (i % 3) * 0.4  # Stack amplicons
            
            # Draw amplicon
            ax.add_patch(Rectangle((amp['start'], y_offset), 
                                   amp['size'], 0.3,
                                   facecolor=color, edgecolor='black', 
                                   linewidth=0.5, alpha=0.7))
            
            # Add label if amplicon is big enough
            if amp['size'] > locus_length * 0.05:
                mid_point = amp['start'] + amp['size'] / 2
                ax.text(mid_point, y_offset + 0.15, 
                       amp['gene'].split('_')[0] if '_' in amp['gene'] else amp['gene'],
                       ha='center', va='center', fontsize=6, rotation=0)
        
        # Calculate and show coverage
        covered_positions = set()
        for amp in amplicons:
            covered_positions.update(range(amp['start'], amp['end']))
        
        coverage_pct = (len(covered_positions) / locus_length) * 100
        
        # Mark uncovered regions in red
        all_positions = set(range(locus_length))
        uncovered = sorted(all_positions - covered_positions)
        
        # Group consecutive uncovered positions into gaps
        if uncovered:
            gaps = []
            gap_start = uncovered[0]
            prev_pos = uncovered[0]
            
            for pos in uncovered[1:]:
                if pos != prev_pos + 1:
                    gaps.append((gap_start, prev_pos + 1))
                    gap_start = pos
                prev_pos = pos
            gaps.append((gap_start, prev_pos + 1))
            
            # Draw gaps
            for gap_start, gap_end in gaps:
                gap_size = gap_end - gap_start
                if gap_size > 100:  # Only show significant gaps
                    ax.add_patch(Rectangle((gap_start, 0), gap_size, 0.5,
                                           facecolor='red', edgecolor='darkred', 
                                           linewidth=1, alpha=0.5))
        
        ax.set_xlim(0, locus_length)
        ax.set_ylim(-0.2, 2.5)
        ax.set_xlabel('Position (bp)', fontsize=12, fontweight='bold')
        ax.set_title(f'{serotype} CPS Locus Coverage\n'
                    f'{len(amplicons)} amplicons, {coverage_pct:.1f}% coverage, '
                    f'{locus_length:,} bp total',
                    fontsize=13, fontweight='bold', pad=15)
        
        # Format x-axis
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))
        
        # Remove y-axis
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        # Add legend
        legend_elements = [
            patches.Patch(facecolor='lightgray', edgecolor='black', label='CPS Locus'),
            patches.Patch(facecolor='blue', alpha=0.7, edgecolor='black', label='Amplicon'),
            patches.Patch(facecolor='red', alpha=0.5, edgecolor='darkred', label='Gap (>100bp)')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
        
        plt.tight_layout()
        
        output_file = self.output_dir / "plots" / f"coverage_{serotype}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved coverage plot: {output_file}")
        
        return output_file
    
    def plot_all_loci_summary(self, max_serotypes=30):
        """Create summary plot showing all loci coverage"""
        logger.info("Creating summary coverage plot...")
        
        # Sort by coverage
        sorted_summary = self.coverage_summary.sort_values('coverage_percent', ascending=False)
        
        # Take top N serotypes
        plot_data = sorted_summary.head(max_serotypes)
        
        fig, ax = plt.subplots(figsize=(14, 10))
        
        y_positions = range(len(plot_data))
        
        for i, (_, row) in enumerate(plot_data.iterrows()):
            serotype = row['serotype']
            locus_length = row['locus_length']
            coverage_pct = row['coverage_percent']
            
            # Draw full locus bar (gray)
            ax.barh(i, locus_length, height=0.8, 
                   color='lightgray', edgecolor='black', linewidth=0.5)
            
            # Draw covered portion (green)
            covered_length = locus_length * (coverage_pct / 100)
            ax.barh(i, covered_length, height=0.8,
                   color='#2E7D32', edgecolor='black', linewidth=0.5, alpha=0.8)
            
            # Add percentage label
            ax.text(locus_length + 500, i, f'{coverage_pct:.1f}%', 
                   va='center', fontsize=9, fontweight='bold')
        
        ax.set_yticks(y_positions)
        ax.set_yticklabels(plot_data['serotype'], fontsize=9)
        ax.set_xlabel('Locus Length (bp)', fontsize=12, fontweight='bold')
        ax.set_title('CPS Locus Coverage by Serotype\n(Top 30)', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Format x-axis
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x/1000):.0f}k'))
        
        # Add legend
        legend_elements = [
            patches.Patch(facecolor='#2E7D32', alpha=0.8, label='Covered by amplicons'),
            patches.Patch(facecolor='lightgray', label='Not covered')
        ]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
        
        plt.tight_layout()
        
        output_file = self.output_dir / "plots" / "coverage_summary_all.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved summary plot: {output_file}")
        
        return output_file
    
    def export_coverage_report(self):
        """Export detailed coverage report"""
        logger.info("Exporting coverage reports...")
        
        # Summary statistics
        summary_file = self.output_dir / "coverage_reports" / "coverage_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("CPS LOCUS COVERAGE SUMMARY\n")
            f.write("=" * 80 + "\n\n")
            
            avg_coverage = self.coverage_summary['coverage_percent'].mean()
            min_coverage = self.coverage_summary['coverage_percent'].min()
            max_coverage = self.coverage_summary['coverage_percent'].max()
            
            f.write(f"Total serotypes: {len(self.coverage_summary)}\n")
            f.write(f"Average coverage: {avg_coverage:.1f}%\n")
            f.write(f"Min coverage: {min_coverage:.1f}%\n")
            f.write(f"Max coverage: {max_coverage:.1f}%\n\n")
            
            # Coverage bins
            bins = [(0, 25), (25, 50), (50, 75), (75, 90), (90, 100)]
            f.write("Coverage distribution:\n")
            for low, high in bins:
                count = len(self.coverage_summary[
                    (self.coverage_summary['coverage_percent'] >= low) & 
                    (self.coverage_summary['coverage_percent'] < high)
                ])
                f.write(f"  {low:3d}-{high:3d}%: {count:3d} serotypes\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write("PER-SEROTYPE DETAILS\n")
            f.write("=" * 80 + "\n\n")
            
            for _, row in self.coverage_summary.sort_values('coverage_percent', ascending=False).iterrows():
                serotype = row['serotype']
                f.write(f"{serotype:20s} | {row['locus_length']:6,d} bp | "
                       f"{row['num_amplicons']:3d} amplicons | "
                       f"{row['coverage_percent']:6.1f}% covered\n")
        
        logger.info(f"Saved summary report: {summary_file}")
        
        # Detailed amplicon mapping
        detail_file = self.output_dir / "coverage_reports" / "amplicon_mapping.csv"
        
        detail_data = []
        for serotype, amplicons in self.amplicon_coverage.items():
            for amp in amplicons:
                detail_data.append({
                    'serotype': serotype,
                    'primer_name': amp['primer_name'],
                    'gene': amp['gene'],
                    'start': amp['start'],
                    'end': amp['end'],
                    'size': amp['size'],
                    'mismatches': amp['total_mismatches']
                })
        
        detail_df = pd.DataFrame(detail_data)
        detail_df.to_csv(detail_file, index=False)
        
        logger.info(f"Saved detailed mapping: {detail_file}")
        
        # Save coverage summary table
        summary_csv = self.output_dir / "coverage_reports" / "coverage_summary.csv"
        self.coverage_summary.to_csv(summary_csv, index=False)
        logger.info(f"Saved coverage summary table: {summary_csv}")
        
        return summary_file


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Visualize amplicon coverage across CPS loci',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze coverage for all serotypes
  python visualize_cps_coverage.py --database ./database --primers primers.csv --output ./coverage_viz
  
  # Plot specific serotypes
  python visualize_cps_coverage.py --database ./database --primers primers.csv --output ./coverage_viz --plot-serotypes 19F,19A,6B
  
  # Allow more mismatches
  python visualize_cps_coverage.py --database ./database --primers primers.csv --output ./coverage_viz --max-mismatches 3
        """
    )
    
    parser.add_argument('--database', required=True, help='Database directory with reference.fasta files')
    parser.add_argument('--primers', required=True, help='Primer CSV/TSV file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--max-mismatches', type=int, default=2, 
                       help='Maximum mismatches allowed for primer binding (default: 2)')
    parser.add_argument('--plot-serotypes', type=str, default=None,
                       help='Comma-separated list of serotypes to plot individually (e.g., "19F,19A,6B")')
    parser.add_argument('--plot-all', action='store_true',
                       help='Plot individual coverage for all serotypes (slow!)')
    parser.add_argument('--debug', action='store_true')
    
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    visualizer = CPSCoverageVisualizer(
        database_dir=args.database,
        primer_file=args.primers,
        output_dir=args.output,
        max_mismatches=args.max_mismatches
    )
    
    try:
        logger.info("=" * 80)
        logger.info("CPS LOCUS COVERAGE ANALYSIS")
        logger.info("=" * 80)
        
        # Load data
        visualizer.load_primers()
        visualizer.load_reference_loci()
        
        # Map amplicons
        visualizer.map_amplicons_to_loci()
        
        # Create summary plots
        visualizer.plot_all_loci_summary()
        
        # Export reports
        visualizer.export_coverage_report()
        
        # Plot individual serotypes if requested
        if args.plot_serotypes:
            serotypes = [s.strip() for s in args.plot_serotypes.split(',')]
            logger.info(f"\nPlotting individual coverage for: {', '.join(serotypes)}")
            for serotype in serotypes:
                visualizer.plot_locus_coverage(serotype)
        
        elif args.plot_all:
            logger.info("\nPlotting all serotypes...")
            for serotype in visualizer.reference_loci.keys():
                visualizer.plot_locus_coverage(serotype)
        
        logger.info("=" * 80)
        logger.info("SUCCESS!")
        logger.info(f"Results saved to: {args.output}")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Failed: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
