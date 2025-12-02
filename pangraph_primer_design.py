#!/usr/bin/env python3
"""
Pan-Genome Graph Primer Designer
Uses k-mer clustering to identify conserved regions across serotypes
and designs universal primers for efficient coverage
"""

import os
import sys
from pathlib import Path
from collections import defaultdict, Counter
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import subprocess
import primer3
import pandas as pd
import numpy as np
import logging
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import itertools

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class PanGenomeGraphPrimerDesigner:
    """Design universal primers using k-mer graph clustering"""
    
    def __init__(self, database_dir, output_dir, 
                 window_size=1000, step_size=500,
                 kmer_size=31, similarity_threshold=0.80,
                 amplicon_size=500, overlap=75):
        self.database_dir = Path(database_dir)
        self.output_dir = Path(output_dir)
        self.window_size = window_size
        self.step_size = step_size
        self.kmer_size = kmer_size
        self.similarity_threshold = similarity_threshold
        self.amplicon_size = amplicon_size
        self.overlap = overlap
        self.min_amplicon = amplicon_size - 100
        self.max_amplicon = amplicon_size + 100
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "primers").mkdir(exist_ok=True)
        (self.output_dir / "clusters").mkdir(exist_ok=True)
        (self.output_dir / "alignments").mkdir(exist_ok=True)
        
        self.reference_loci = {}
        self.windows = []
        self.clusters = []
        self.primers = []
        self.primer_pools = {'pool1': [], 'pool2': []}
        
    def load_reference_loci(self):
        """Load all reference CPS loci"""
        logger.info("Loading reference CPS loci...")
        
        # Root reference.fasta
        root_ref = self.database_dir / "reference.fasta"
        if root_ref.exists():
            for record in SeqIO.parse(root_ref, "fasta"):
                self.reference_loci[record.id] = {
                    'id': record.id,
                    'sequence': str(record.seq).upper(),
                    'length': len(record.seq)
                }
        
        # Subdirectory reference.fasta files
        for subdir in self.database_dir.iterdir():
            if subdir.is_dir():
                ref_file = subdir / "reference.fasta"
                if ref_file.exists():
                    for record in SeqIO.parse(ref_file, "fasta"):
                        serotype = subdir.name
                        if serotype not in self.reference_loci:
                            self.reference_loci[serotype] = {
                                'id': record.id,
                                'sequence': str(record.seq).upper(),
                                'length': len(record.seq)
                            }
        
        logger.info(f"Loaded {len(self.reference_loci)} reference loci")
        
        total_bp = sum(data['length'] for data in self.reference_loci.values())
        logger.info(f"Total sequence: {total_bp:,} bp")
        
        return self.reference_loci
    
    def extract_windows(self):
        """Extract overlapping windows from all loci"""
        logger.info(f"\nExtracting {self.window_size}bp windows (step={self.step_size}bp)...")
        
        window_id = 0
        
        for serotype, locus_data in self.reference_loci.items():
            sequence = locus_data['sequence']
            seq_len = len(sequence)
            
            num_windows = 0
            for start in range(0, seq_len - self.window_size + 1, self.step_size):
                end = start + self.window_size
                window_seq = sequence[start:end]
                
                # Skip windows with too many N's
                n_count = window_seq.count('N')
                if n_count > self.window_size * 0.1:  # >10% N's
                    continue
                
                self.windows.append({
                    'window_id': f'win_{window_id}',
                    'serotype': serotype,
                    'start': start,
                    'end': end,
                    'sequence': window_seq,
                    'length': len(window_seq)
                })
                
                window_id += 1
                num_windows += 1
            
            logger.debug(f"  {serotype}: {num_windows} windows")
        
        logger.info(f"Extracted {len(self.windows)} windows total")
        
        return self.windows
    
    def _get_kmer_profile(self, sequence):
        """Get k-mer profile"""
        kmers = [sequence[i:i+self.kmer_size] 
                for i in range(len(sequence) - self.kmer_size + 1)]
        # Filter k-mers with N's
        kmers = [km for km in kmers if 'N' not in km]
        return Counter(kmers)
    
    def _kmer_jaccard_similarity(self, profile1, profile2):
        """Calculate Jaccard similarity"""
        if not profile1 or not profile2:
            return 0.0
        
        kmers1 = set(profile1.keys())
        kmers2 = set(profile2.keys())
        
        intersection = len(kmers1 & kmers2)
        union = len(kmers1 | kmers2)
        
        return intersection / union if union > 0 else 0.0
    
    def cluster_windows(self, min_cluster_size=1):
        """Cluster similar windows using k-mer profiles"""
        logger.info(f"\nClustering windows (k={self.kmer_size}, threshold={self.similarity_threshold})...")
        logger.info(f"  Min cluster size: {min_cluster_size} (including unique regions)")
        
        # Compute k-mer profiles for all windows
        logger.info("Computing k-mer profiles...")
        for window in self.windows:
            window['kmer_profile'] = self._get_kmer_profile(window['sequence'])
        
        # Build similarity matrix (sample if too large)
        n = len(self.windows)
        
        if n > 5000:
            logger.warning(f"Large dataset ({n} windows). Using sparse clustering...")
            # Use all windows but sparse computation
            working_windows = self.windows
            use_sparse = True
        else:
            working_windows = self.windows
            use_sparse = False
        
        logger.info("Building similarity matrix...")
        n_work = len(working_windows)
        
        # Use sparse approach - only compute for likely similar pairs
        similarities = {}
        
        for i in range(n_work):
            if (i + 1) % 500 == 0:
                logger.info(f"  Processing window {i+1}/{n_work}...")
            
            for j in range(i + 1, min(i + 100, n_work)):  # Only compare nearby windows
                sim = self._kmer_jaccard_similarity(
                    working_windows[i]['kmer_profile'],
                    working_windows[j]['kmer_profile']
                )
                
                # Only store if above threshold
                if sim >= self.similarity_threshold:
                    similarities[(i, j)] = sim
        
        logger.info(f"Found {len(similarities)} similar pairs")
        
        # Build clusters using connected components
        logger.info("Building clusters from similar pairs...")
        clusters_dict = {}
        cluster_id = 0
        
        for i in range(n_work):
            clusters_dict[i] = cluster_id
            cluster_id += 1
        
        # Merge clusters based on similarities
        for (i, j), sim in similarities.items():
            # Merge clusters
            cluster_i = clusters_dict[i]
            cluster_j = clusters_dict[j]
            
            if cluster_i != cluster_j:
                # Merge j's cluster into i's cluster
                old_cluster = cluster_j
                for idx in range(n_work):
                    if clusters_dict[idx] == old_cluster:
                        clusters_dict[idx] = cluster_i
        
        # Group windows by cluster
        cluster_groups = defaultdict(list)
        for idx, cluster in clusters_dict.items():
            cluster_groups[cluster].append(working_windows[idx])
        
        # Filter by minimum size (but keep singletons if min_cluster_size=1)
        self.clusters = []
        singletons = 0
        for cluster_id, windows in cluster_groups.items():
            if len(windows) >= min_cluster_size:
                self.clusters.append({
                    'cluster_id': f'cluster_{len(self.clusters)}',
                    'windows': windows,
                    'num_serotypes': len(set(w['serotype'] for w in windows))
                })
                if len(windows) == 1:
                    singletons += 1
        
        logger.info(f"\nCreated {len(self.clusters)} clusters ({singletons} unique regions)")
        
        # Show top clusters
        sorted_clusters = sorted(self.clusters, key=lambda x: len(x['windows']), reverse=True)
        logger.info("\nTop universal clusters:")
        for i, cluster in enumerate(sorted_clusters[:20], 1):
            serotypes = set(w['serotype'] for w in cluster['windows'])
            serotype_str = ', '.join(list(serotypes)[:3])
            if len(serotypes) > 3:
                serotype_str += f' ...+{len(serotypes)-3}'
            
            logger.info(f"  [{i:2d}] {cluster['cluster_id']}: "
                       f"{len(cluster['windows'])} windows from {len(serotypes)} serotypes | {serotype_str}")
        
        return self.clusters
    
    def design_primers_for_clusters(self):
        """Design primers for each cluster"""
        logger.info("\nDesigning primers for clusters...")
        
        for i, cluster in enumerate(self.clusters, 1):
            cluster_id = cluster['cluster_id']
            windows = cluster['windows']
            
            logger.info(f"[{i}/{len(self.clusters)}] {cluster_id} "
                       f"({len(windows)} windows, {cluster['num_serotypes']} serotypes)")
            
            # Get representative sequence (longest or consensus)
            if len(windows) == 1:
                consensus = windows[0]['sequence']
            elif len(windows) <= 10:
                # Align and make consensus for small clusters
                consensus = self._make_consensus_from_windows(windows)
            else:
                # Use longest for large clusters
                consensus = max(windows, key=lambda x: x['length'])['sequence']
            
            # Design tiling primers across this region
            primers = self._design_tiled_primers(
                consensus, 
                cluster_id, 
                len(windows),
                cluster['num_serotypes']
            )
            
            if primers:
                for primer in primers:
                    primer['cluster_id'] = cluster_id
                    primer['num_windows'] = len(windows)
                    primer['num_serotypes'] = cluster['num_serotypes']
                    primer['serotypes_covered'] = ','.join(sorted(set(w['serotype'] for w in windows)))
                    primer['primer_type'] = 'universal'
                
                self.primers.extend(primers)
                logger.info(f"  ✓ {len(primers)} primers")
            else:
                logger.warning(f"  ✗ No primers")
        
        logger.info(f"\nTotal universal primers designed: {len(self.primers)}")
        
        return self.primers
    
    def fill_coverage_gaps(self, max_mismatches=2):
        """Design primers to fill gaps in coverage"""
        logger.info("\n" + "=" * 80)
        logger.info("FILLING COVERAGE GAPS")
        logger.info("=" * 80)
        
        gap_primers = []
        
        # First, map which windows are covered by existing clusters
        window_coverage = defaultdict(set)  # serotype -> set of (start, end) tuples
        
        for cluster in self.clusters:
            for window in cluster['windows']:
                serotype = window['serotype']
                window_coverage[serotype].add((window['start'], window['end']))
        
        logger.info(f"Mapped coverage from {len(self.clusters)} clusters")
        
        # For each serotype, find gaps
        for serotype, locus_data in sorted(self.reference_loci.items()):
            sequence = locus_data['sequence']
            seq_len = len(sequence)
            
            # Mark covered positions
            covered = [False] * seq_len
            
            if serotype in window_coverage:
                for start, end in window_coverage[serotype]:
                    for i in range(start, min(end, seq_len)):
                        covered[i] = True
            
            # Find gaps
            gaps = []
            in_gap = False
            gap_start = 0
            
            for i in range(seq_len):
                if not covered[i]:
                    if not in_gap:
                        gap_start = i
                        in_gap = True
                else:
                    if in_gap:
                        gap_len = i - gap_start
                        if gap_len >= self.min_amplicon:  # Only fill significant gaps
                            gaps.append((gap_start, i))
                        in_gap = False
            
            if in_gap and (seq_len - gap_start) >= self.min_amplicon:
                gaps.append((gap_start, seq_len))
            
            # Calculate coverage
            covered_bp = sum(covered)
            coverage_pct = (covered_bp / seq_len) * 100
            
            logger.info(f"\n{serotype}: {coverage_pct:.1f}% covered by universal primers")
            
            if gaps:
                logger.info(f"  Found {len(gaps)} gaps totaling {sum(e-s for s,e in gaps):,} bp")
                
                # Design primers for gaps
                gap_primers_for_serotype = []
                amp_num = 0
                
                for gap_start, gap_end in gaps:
                    gap_len = gap_end - gap_start
                    gap_seq = sequence[gap_start:gap_end]
                    
                    logger.info(f"  Gap: {gap_start:,}-{gap_end:,} ({gap_len:,} bp)")
                    
                    # Tile this gap
                    step = self.amplicon_size - self.overlap
                    
                    for region_start in range(0, gap_len, step):
                        region_end = min(region_start + self.amplicon_size + 200, gap_len)
                        
                        if region_end - region_start < self.min_amplicon:
                            break
                        
                        amp_num += 1
                        
                        target_seq = gap_seq[region_start:region_end]
                        start = region_start + 100
                        end = min(start + self.amplicon_size, region_end)
                        
                        if len(target_seq) < 100:
                            continue
                        
                        clean_seq = ''.join([b if b in 'ACGTN' else 'N' for b in target_seq])
                        
                        result = self._design_pair(
                            clean_seq, amp_num, f"{serotype}_gap", 1,
                            start - region_start, end - region_start, 
                            gap_start + region_start
                        )
                        
                        if result:
                            result['serotype'] = serotype
                            result['cluster_id'] = f'{serotype}_gap'
                            result['num_serotypes'] = 1
                            result['serotypes_covered'] = serotype
                            result['primer_type'] = 'gap_fill'
                            result['gap_region'] = f'{gap_start}-{gap_end}'
                            gap_primers_for_serotype.append(result)
                
                if gap_primers_for_serotype:
                    logger.info(f"  Designed {len(gap_primers_for_serotype)} gap-filling primers")
                    gap_primers.extend(gap_primers_for_serotype)
            else:
                logger.info(f"  No gaps - fully covered!")
        
        logger.info(f"\n{'='*80}")
        logger.info(f"Total gap-filling primers: {len(gap_primers)}")
        logger.info(f"Total primers (universal + gap-fill): {len(self.primers) + len(gap_primers)}")
        
        # Report coverage improvement
        total_locus_bp = sum(data['length'] for data in self.reference_loci.values())
        naive_primers = total_locus_bp / self.amplicon_size
        total_primers = len(self.primers) + len(gap_primers)
        reduction = ((naive_primers - total_primers) / naive_primers) * 100
        
        logger.info(f"Efficiency: {reduction:.1f}% reduction vs naive tiling")
        logger.info(f"  Naive approach: ~{int(naive_primers)} primers")
        logger.info(f"  Pan-graph + gap-fill: {total_primers} primers")
        logger.info(f"{'='*80}")
        
        self.primers.extend(gap_primers)
        
        return gap_primers
    
    def _make_consensus_from_windows(self, windows):
        """Make consensus from multiple window sequences"""
        if len(windows) == 1:
            return windows[0]['sequence']
        
        # Write to temp file
        temp_fasta = self.output_dir / "temp_align.fasta"
        with open(temp_fasta, 'w') as f:
            for i, win in enumerate(windows):
                f.write(f">{i}\n{win['sequence']}\n")
        
        try:
            # Align with MAFFT
            cmd = ['mafft', '--quiet', '--auto', str(temp_fasta)]
            result = subprocess.run(cmd, capture_output=True, text=True, 
                                  check=True, timeout=60)
            
            # Parse alignment
            alignment = AlignIO.read(result.stdout.split('\n'), "fasta")
            
            # Make consensus
            consensus = []
            for i in range(alignment.get_alignment_length()):
                column = alignment[:, i]
                bases = [b.upper() for b in column if b.upper() in 'ACGT']
                
                if len(bases) > 0:
                    most_common = Counter(bases).most_common(1)[0]
                    if most_common[1] / len(bases) > 0.5:
                        consensus.append(most_common[0])
                    else:
                        consensus.append('N')
            
            return ''.join(consensus)
        
        except Exception as e:
            logger.debug(f"Alignment failed: {e}")
            # Fallback to longest sequence
            return max(windows, key=lambda x: x['length'])['sequence']
        
        finally:
            if temp_fasta.exists():
                temp_fasta.unlink()
    
    def _design_tiled_primers(self, sequence, cluster_id, num_windows, num_serotypes):
        """Design tiling primers across a sequence"""
        primers = []
        seq_len = len(sequence)
        step = self.amplicon_size - self.overlap
        
        amp_num = 0
        for region_start in range(0, seq_len, step):
            region_end = min(region_start + self.amplicon_size + 200, seq_len)
            
            if region_end - region_start < self.min_amplicon:
                break
            
            amp_num += 1
            
            target_seq = sequence[region_start:region_end]
            start = region_start + 100
            end = min(start + self.amplicon_size, region_end)
            
            if len(target_seq) < 100:
                continue
            
            clean_seq = ''.join([b if b in 'ACGTN' else 'N' for b in target_seq])
            
            result = self._design_pair(
                clean_seq, amp_num, cluster_id, num_windows,
                start - region_start, end - region_start, region_start
            )
            
            if result:
                primers.append(result)
        
        return primers
    
    def _design_pair(self, sequence, amplicon_num, gene_name, num_alleles,
                    target_start, target_end, offset):
        """Design single primer pair"""
        
        seq_args = {
            'SEQUENCE_ID': f'{gene_name}_amp{amplicon_num}',
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_TARGET': [target_start, target_end - target_start],
        }
        
        global_args = {
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_OPT_SIZE': 25,
            'PRIMER_MIN_SIZE': 22,
            'PRIMER_MAX_SIZE': 30,
            'PRIMER_OPT_TM': 61.0,
            'PRIMER_MIN_TM': 59.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[self.min_amplicon, self.max_amplicon]],
            'PRIMER_MAX_POLY_X': 4,
            'PRIMER_MAX_NS_ACCEPTED': 1,
            'PRIMER_MAX_SELF_ANY': 8,
            'PRIMER_MAX_SELF_END': 3,
            'PRIMER_PAIR_MAX_COMPL_ANY': 8,
            'PRIMER_PAIR_MAX_COMPL_END': 3,
            'PRIMER_NUM_RETURN': 5,
        }
        
        try:
            result = primer3.bindings.design_primers(seq_args, global_args)
            
            if result.get('PRIMER_PAIR_NUM_RETURNED', 0) > 0:
                return {
                    'gene': gene_name,
                    'amplicon_num': amplicon_num,
                    'num_alleles': num_alleles,
                    'name': f'{gene_name}_amp{amplicon_num}',
                    'left_seq': result['PRIMER_LEFT_0_SEQUENCE'],
                    'right_seq': result['PRIMER_RIGHT_0_SEQUENCE'],
                    'left_tm': result['PRIMER_LEFT_0_TM'],
                    'right_tm': result['PRIMER_RIGHT_0_TM'],
                    'left_gc': result['PRIMER_LEFT_0_GC_PERCENT'],
                    'right_gc': result['PRIMER_RIGHT_0_GC_PERCENT'],
                    'product_size': result['PRIMER_PAIR_0_PRODUCT_SIZE'],
                    'left_pos': result['PRIMER_LEFT_0'][0] + offset,
                    'right_pos': result['PRIMER_RIGHT_0'][0] + offset,
                    'penalty': result.get('PRIMER_PAIR_0_PENALTY', 0)
                }
        except Exception as e:
            logger.debug(f"Primer design exception: {e}")
        
        return None
    
    def assign_primer_pools(self):
        """Assign primers to pools"""
        logger.info("\nAssigning primers to pools...")
        
        for primer in self.primers:
            if primer['amplicon_num'] % 2 == 0:
                self.primer_pools['pool1'].append(primer)
                primer['pool'] = 'pool1'
            else:
                self.primer_pools['pool2'].append(primer)
                primer['pool'] = 'pool2'
        
        logger.info(f"  Pool 1: {len(self.primer_pools['pool1'])} primers")
        logger.info(f"  Pool 2: {len(self.primer_pools['pool2'])} primers")
        
        return self.primer_pools
    
    def export_primers(self, format='csv'):
        """Export primers"""
        logger.info("\nExporting primers...")
        
        if not self.primers:
            logger.error("No primers to export!")
            return None
        
        df = pd.DataFrame(self.primers)
        df['left_name'] = df.apply(lambda x: f"{x['name']}_LEFT", axis=1)
        df['right_name'] = df.apply(lambda x: f"{x['name']}_RIGHT", axis=1)
        
        cols = ['gene', 'cluster_id', 'amplicon_num', 'num_alleles', 'num_serotypes',
                'serotypes_covered', 'pool', 'name',
                'left_name', 'left_seq', 'left_tm', 'left_gc', 'left_pos',
                'right_name', 'right_seq', 'right_tm', 'right_gc', 'right_pos',
                'product_size', 'penalty']
        
        df = df[[col for col in cols if col in df.columns]]
        
        sep = '\t' if format == 'tsv' else ','
        out = self.output_dir / "primers" / f"pangraph_primers.{format}"
        df.to_csv(out, index=False, sep=sep)
        
        logger.info(f"  All primers: {out}")
        
        # Export pools
        for pool_name, pool_primers in self.primer_pools.items():
            if pool_primers:
                pool_df = pd.DataFrame(pool_primers)
                pool_df['left_name'] = pool_df.apply(lambda x: f"{x['name']}_LEFT", axis=1)
                pool_df['right_name'] = pool_df.apply(lambda x: f"{x['name']}_RIGHT", axis=1)
                
                pool_df = pool_df[[col for col in cols if col in pool_df.columns]]
                pool_file = self.output_dir / "primers" / f"{pool_name}_primers.{format}"
                pool_df.to_csv(pool_file, index=False, sep=sep)
                logger.info(f"  {pool_name}: {pool_file}")
        
        return out
    
    def generate_report(self):
        """Generate report"""
        report_file = self.output_dir / "pangraph_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("PAN-GENOME GRAPH PRIMER DESIGN REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("PARAMETERS:\n")
            f.write(f"  Window size: {self.window_size} bp\n")
            f.write(f"  Step size: {self.step_size} bp\n")
            f.write(f"  K-mer size: {self.kmer_size}\n")
            f.write(f"  Similarity threshold: {self.similarity_threshold}\n")
            f.write(f"  Amplicon size: {self.amplicon_size} bp\n\n")
            
            f.write("RESULTS:\n")
            f.write(f"  Reference loci: {len(self.reference_loci)}\n")
            f.write(f"  Windows extracted: {len(self.windows)}\n")
            f.write(f"  Clusters formed: {len(self.clusters)}\n")
            f.write(f"  Universal primers designed: {len(self.primers)}\n")
            f.write(f"  Pool 1: {len(self.primer_pools['pool1'])}\n")
            f.write(f"  Pool 2: {len(self.primer_pools['pool2'])}\n\n")
            
            # Efficiency metrics
            total_locus_bp = sum(data['length'] for data in self.reference_loci.values())
            naive_primers = total_locus_bp / self.amplicon_size
            reduction = ((naive_primers - len(self.primers)) / naive_primers) * 100
            
            f.write("EFFICIENCY:\n")
            f.write(f"  Total sequence: {total_locus_bp:,} bp\n")
            f.write(f"  Naive approach: ~{int(naive_primers)} primers\n")
            f.write(f"  Pan-graph approach: {len(self.primers)} primers\n")
            f.write(f"  Reduction: {reduction:.1f}%\n\n")
            
            f.write("TOP UNIVERSAL PRIMERS (covering most serotypes):\n")
            sorted_primers = sorted(self.primers, 
                                   key=lambda x: x.get('num_serotypes', 0), 
                                   reverse=True)
            for i, primer in enumerate(sorted_primers[:20], 1):
                f.write(f"  {i:2d}. {primer['name']:30s} "
                       f"covers {primer.get('num_serotypes', 0):2d} serotypes\n")
            
            f.write("\n" + "=" * 80 + "\n")
        
        logger.info(f"\nReport saved: {report_file}")
        return report_file


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Pan-genome graph-based universal primer design',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--database', required=True, help='Database directory')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--window-size', type=int, default=1000,
                       help='Window size for segmentation (bp, default: 1000)')
    parser.add_argument('--step-size', type=int, default=500,
                       help='Step size between windows (bp, default: 500)')
    parser.add_argument('--kmer-size', type=int, default=31,
                       help='K-mer size for clustering (default: 31)')
    parser.add_argument('--similarity', type=float, default=0.80,
                       help='Similarity threshold for clustering (default: 0.80)')
    parser.add_argument('--amplicon-size', type=int, default=500,
                       help='Target amplicon size (bp, default: 500)')
    parser.add_argument('--overlap', type=int, default=75,
                       help='Amplicon overlap (bp, default: 75)')
    parser.add_argument('--format', choices=['csv', 'tsv'], default='csv')
    parser.add_argument('--debug', action='store_true')
    
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    designer = PanGenomeGraphPrimerDesigner(
        database_dir=args.database,
        output_dir=args.output,
        window_size=args.window_size,
        step_size=args.step_size,
        kmer_size=args.kmer_size,
        similarity_threshold=args.similarity,
        amplicon_size=args.amplicon_size,
        overlap=args.overlap
    )
    
    try:
        logger.info("=" * 80)
        logger.info("PAN-GENOME GRAPH PRIMER DESIGN")
        logger.info("=" * 80)
        
        # Load and process
        designer.load_reference_loci()
        designer.extract_windows()
        designer.cluster_windows(min_cluster_size=1)  # Include unique regions
        designer.design_primers_for_clusters()
        
        # Fill any coverage gaps with serotype-specific primers
        logger.info("\n" + "=" * 80)
        logger.info("PHASE 2: GAP FILLING")
        logger.info("=" * 80)
        designer.fill_coverage_gaps()
        
        if not designer.primers:
            logger.error("No primers designed!")
            sys.exit(1)
        
        designer.assign_primer_pools()
        designer.export_primers(format=args.format)
        designer.generate_report()
        
        logger.info("=" * 80)
        logger.info("SUCCESS!")
        logger.info(f"Designed {len(designer.primers)} universal primers")
        logger.info(f"Results: {args.output}")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Failed: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()