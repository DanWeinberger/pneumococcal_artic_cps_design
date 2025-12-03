#!/usr/bin/env python3
"""
In Silico Mixture Validation
Simulates mixed serotype samples, generates reads, and evaluates serotype calling accuracy
"""

import os
import sys
from pathlib import Path
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import logging
import random
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class MixtureSimulator:
    """Simulate mixed serotype samples and evaluate primer performance"""
    
    def __init__(self, database_dir, primer_file, output_dir, 
                 max_mismatches=2, read_length=150, coverage_depth=100):
        self.database_dir = Path(database_dir)
        self.primer_file = Path(primer_file)
        self.output_dir = Path(output_dir)
        self.max_mismatches = max_mismatches
        self.read_length = read_length
        self.coverage_depth = coverage_depth
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "simulations").mkdir(exist_ok=True)
        (self.output_dir / "results").mkdir(exist_ok=True)
        (self.output_dir / "plots").mkdir(exist_ok=True)
        
        self.reference_loci = {}
        self.primers = None
        self.amplicon_map = {}
        
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
    
    def map_primers(self):
        """Map primers to all serotypes"""
        logger.info("\nMapping primers to reference loci...")
        
        for primer_row in self.primers.itertuples():
            primer_name = primer_row.name
            left_seq = primer_row.left_seq
            right_seq = primer_row.right_seq
            right_seq_rc = str(Seq(right_seq).reverse_complement())
            
            self.amplicon_map[primer_name] = {}
            
            for serotype, sequence in self.reference_loci.items():
                left_match = self._find_primer_in_sequence(left_seq, sequence, self.max_mismatches)
                right_match = self._find_primer_in_sequence(right_seq_rc, sequence, self.max_mismatches)
                
                if left_match and right_match and right_match['position'] > left_match['position']:
                    amplicon_start = left_match['position']
                    amplicon_end = right_match['position'] + len(right_seq)
                    amplicon_size = amplicon_end - amplicon_start
                    
                    if 200 <= amplicon_size <= 1000:
                        self.amplicon_map[primer_name][serotype] = {
                            'start': amplicon_start,
                            'end': amplicon_end,
                            'amplicon_seq': sequence[amplicon_start:amplicon_end]
                        }
        
        logger.info(f"Mapped {len(self.amplicon_map)} primers")
        return self.amplicon_map
    
    def simulate_pcr_amplification(self, serotypes, proportions, num_pcr_cycles=25):
        """
        Simulate PCR amplification with exponential growth
        
        Args:
            serotypes: List of serotype names in mixture
            proportions: Initial proportions of each serotype
            num_pcr_cycles: Number of PCR cycles (default: 25)
        
        Returns:
            Dictionary of amplicon counts after PCR
        """
        logger.info(f"\nSimulating PCR amplification ({num_pcr_cycles} cycles)...")
        
        # Initial template molecules (before PCR)
        initial_molecules = 1000  # Starting template molecules
        template_counts = {}
        
        for serotype, proportion in zip(serotypes, proportions):
            template_counts[serotype] = int(initial_molecules * proportion)
        
        logger.info(f"  Initial template molecules:")
        for st, count in template_counts.items():
            logger.info(f"    {st}: {count} molecules")
        
        # Track amplicon counts per primer per serotype
        amplicon_counts = defaultdict(lambda: defaultdict(int))
        
        # For each primer, simulate PCR
        for primer_name in self.amplicon_map.keys():
            # Check which serotypes have this amplicon
            available_serotypes = [st for st in serotypes if st in self.amplicon_map[primer_name]]
            
            if not available_serotypes:
                continue
            
            # Initialize amplicon count = template count
            for serotype in available_serotypes:
                amplicon_counts[primer_name][serotype] = template_counts[serotype]
            
            # Simulate PCR cycles
            for cycle in range(num_pcr_cycles):
                for serotype in available_serotypes:
                    current_count = amplicon_counts[primer_name][serotype]
                    
                    # PCR efficiency (not perfect doubling)
                    # Account for primer binding efficiency and mismatches
                    amplicon_info = self.amplicon_map[primer_name][serotype]
                    
                    # Efficiency decreases with mismatches
                    # Also add some stochasticity
                    base_efficiency = 0.90  # 90% efficiency per cycle
                    efficiency = base_efficiency * np.random.normal(1.0, 0.05)  # Add noise
                    
                    # New molecules this cycle
                    new_molecules = int(current_count * efficiency)
                    amplicon_counts[primer_name][serotype] += new_molecules
        
        # Log final counts
        logger.info(f"\n  Amplicon counts after PCR:")
        total_amplicons = 0
        for primer_name in list(amplicon_counts.keys())[:5]:  # Show first 5
            for serotype, count in amplicon_counts[primer_name].items():
                logger.info(f"    {primer_name} ({serotype}): {count:,} molecules")
                total_amplicons += count
        logger.info(f"  ... (showing first 5 primers)")
        logger.info(f"  Total amplicons: {total_amplicons:,}")
        
        return amplicon_counts
    
    def generate_sequencing_reads(self, amplicon_counts, coverage_depth):
        """
        Generate sequencing reads from PCR products
        Simulates library preparation and sequencing
        
        Args:
            amplicon_counts: Dictionary from simulate_pcr_amplification
            coverage_depth: Target coverage per amplicon
        
        Returns:
            List of reads and metadata
        """
        logger.info(f"\nGenerating sequencing reads (target coverage: {coverage_depth}×)...")
        
        all_reads = []
        read_metadata = []
        
        # Calculate total amplicon pool
        total_amplicons = sum(
            sum(serotype_counts.values()) 
            for serotype_counts in amplicon_counts.values()
        )
        
        if total_amplicons == 0:
            logger.error("No amplicons generated!")
            return [], []
        
        # For each amplicon, generate reads proportional to its abundance
        for primer_name, serotype_counts in amplicon_counts.items():
            for serotype, count in serotype_counts.items():
                # Number of reads proportional to amplicon abundance
                read_fraction = count / total_amplicons
                num_reads = int(coverage_depth * read_fraction * len(amplicon_counts))
                
                if num_reads == 0:
                    continue
                
                # Get amplicon sequence
                amplicon_info = self.amplicon_map[primer_name][serotype]
                amplicon_seq = amplicon_info['amplicon_seq']
                
                # Generate reads from this amplicon
                for _ in range(num_reads):
                    # Random position in amplicon
                    if len(amplicon_seq) <= self.read_length:
                        read_seq = amplicon_seq
                    else:
                        start = random.randint(0, len(amplicon_seq) - self.read_length)
                        read_seq = amplicon_seq[start:start + self.read_length]
                    
                    # Add sequencing errors (Illumina-like error rate)
                    error_rate = 0.001
                    read_list = list(read_seq)
                    for i in range(len(read_list)):
                        if random.random() < error_rate:
                            read_list[i] = random.choice(['A', 'C', 'G', 'T'])
                    
                    final_read = ''.join(read_list)
                    all_reads.append(final_read)
                    
                    read_metadata.append({
                        'read': final_read,
                        'primer': primer_name,
                        'true_serotype': serotype,
                        'amplicon_seq': amplicon_seq,
                        'pcr_count': count
                    })
        
        logger.info(f"  Generated {len(all_reads):,} reads")
        
        return all_reads, read_metadata
    
    def simulate_mixture(self, serotypes, proportions, mixture_name):
        """
        Simulate a mixed sample with PCR amplification and sequencing
        
        Args:
            serotypes: List of serotype names
            proportions: List of proportions (must sum to 1.0)
            mixture_name: Name for this mixture
        """
        logger.info(f"\n{'='*80}")
        logger.info(f"SIMULATING MIXTURE: {mixture_name}")
        logger.info(f"{'='*80}")
        logger.info(f"  Serotypes: {serotypes}")
        logger.info(f"  Proportions: {[f'{p*100:.1f}%' for p in proportions]}")
        
        if abs(sum(proportions) - 1.0) > 0.01:
            raise ValueError("Proportions must sum to 1.0")
        
        # Step 1: Simulate PCR amplification
        amplicon_counts = self.simulate_pcr_amplification(serotypes, proportions)
        
        # Step 2: Generate sequencing reads from PCR products
        reads, read_metadata = self.generate_sequencing_reads(amplicon_counts, self.coverage_depth)
        
        return reads, read_metadata
    
    def build_kmer_index(self, k=31):
        """
        Build k-mer index of all amplicons for fast read mapping
        This is what real aligners (BWA, Bowtie, minimap2) do
        """
        logger.info(f"\nBuilding k-mer index (k={k})...")
        
        kmer_index = defaultdict(list)  # kmer -> list of (primer, serotype, position)
        
        for primer_name, serotype_dict in self.amplicon_map.items():
            for serotype, amplicon_info in serotype_dict.items():
                amplicon_seq = amplicon_info['amplicon_seq']
                
                # Extract all k-mers from this amplicon
                for i in range(len(amplicon_seq) - k + 1):
                    kmer = amplicon_seq[i:i+k]
                    kmer_index[kmer].append({
                        'primer': primer_name,
                        'serotype': serotype,
                        'position': i
                    })
        
        logger.info(f"  Indexed {len(kmer_index):,} unique k-mers")
        return kmer_index
    
    def map_read_with_kmers(self, read_seq, kmer_index, k=31):
        """
        Map a read to amplicons using k-mer seeds (like BWA-MEM)
        
        Strategy:
        1. Extract k-mers from read
        2. Find matching k-mers in index (seeds)
        3. Extend seeds to full alignment
        4. Return best match
        """
        # Extract k-mers from read
        read_kmers = []
        for i in range(len(read_seq) - k + 1):
            kmer = read_seq[i:i+k]
            read_kmers.append((kmer, i))  # (kmer, position in read)
        
        # Find seed matches
        candidate_amplicons = defaultdict(list)  # (primer, serotype) -> list of seed positions
        
        for kmer, read_pos in read_kmers:
            if kmer in kmer_index:
                for hit in kmer_index[kmer]:
                    key = (hit['primer'], hit['serotype'])
                    candidate_amplicons[key].append({
                        'read_pos': read_pos,
                        'amplicon_pos': hit['position']
                    })
        
        if not candidate_amplicons:
            return None  # No k-mer matches found
        
        # Score each candidate amplicon based on number of matching k-mers
        best_match = None
        best_score = 0
        
        for (primer, serotype), seeds in candidate_amplicons.items():
            # Number of matching k-mers = score
            score = len(seeds)
            
            if score > best_score:
                best_score = score
                best_match = (primer, serotype)
        
        return best_match
    
    def map_reads_to_amplicons(self, reads, read_map, use_kmers=True, k=31):
        """
        Map reads back to amplicons to identify serotypes
        
        Args:
            reads: List of read sequences
            read_map: Metadata about reads (for validation)
            use_kmers: If True, use k-mer indexing (FAST). If False, use naive sliding (SLOW)
            k: K-mer size (default: 31)
        """
        logger.info(f"\nMapping reads to amplicons...")
        logger.info(f"  Method: {'K-mer indexing' if use_kmers else 'Naive sliding window'}")
        logger.info(f"  Total reads: {len(reads):,}")
        
        if use_kmers:
            # Build k-mer index once
            kmer_index = self.build_kmer_index(k=k)
        
        # Count reads per primer per serotype
        primer_serotype_counts = defaultdict(lambda: defaultdict(int))
        
        mapped_count = 0
        unmapped_count = 0
        
        for idx, read_info in enumerate(read_map):
            if (idx + 1) % 10000 == 0:
                logger.info(f"    Mapped {idx+1:,}/{len(read_map):,} reads...")
            
            read_seq = read_info['read']
            
            if use_kmers:
                # Fast k-mer based mapping
                match = self.map_read_with_kmers(read_seq, kmer_index, k=k)
                
                if match:
                    primer_name, serotype = match
                    primer_serotype_counts[primer_name][serotype] += 1
                    mapped_count += 1
                else:
                    unmapped_count += 1
            
            else:
                # Slow naive approach (slide across every amplicon)
                best_match = None
                best_mismatches = float('inf')
                
                for primer_name, serotype_dict in self.amplicon_map.items():
                    for serotype, amplicon_info in serotype_dict.items():
                        amplicon_seq = amplicon_info['amplicon_seq']
                        
                        # Check if read matches this amplicon
                        if read_seq in amplicon_seq:
                            mismatches = 0
                        else:
                            # Slide read across amplicon
                            min_mm = float('inf')
                            for i in range(len(amplicon_seq) - len(read_seq) + 1):
                                amplicon_substr = amplicon_seq[i:i+len(read_seq)]
                                mm = sum(1 for a, b in zip(read_seq, amplicon_substr) if a != b)
                                if mm < min_mm:
                                    min_mm = mm
                            mismatches = min_mm
                        
                        if mismatches < best_mismatches:
                            best_mismatches = mismatches
                            best_match = (primer_name, serotype)
                
                if best_match and best_mismatches <= 5:
                    primer_serotype_counts[best_match[0]][best_match[1]] += 1
                    mapped_count += 1
                else:
                    unmapped_count += 1
        
        logger.info(f"\n  Mapped: {mapped_count:,} reads ({mapped_count/len(read_map)*100:.1f}%)")
        logger.info(f"  Unmapped: {unmapped_count:,} reads ({unmapped_count/len(read_map)*100:.1f}%)")
        
        return primer_serotype_counts
    
    def deconvolve_mixture(self, primer_serotype_counts, mixture_serotypes):
        """
        Estimate serotype proportions from read counts
        """
        logger.info("\nDeconvolving mixture...")
        
        # Aggregate across all primers
        serotype_total_reads = defaultdict(int)
        
        for primer_name, serotype_counts in primer_serotype_counts.items():
            for serotype, count in serotype_counts.items():
                serotype_total_reads[serotype] += count
        
        # Calculate proportions
        total_reads = sum(serotype_total_reads.values())
        
        if total_reads == 0:
            logger.error("No reads mapped!")
            return {}
        
        estimated_proportions = {}
        for serotype, count in serotype_total_reads.items():
            estimated_proportions[serotype] = count / total_reads
        
        # Sort by proportion
        estimated_proportions = dict(sorted(estimated_proportions.items(), 
                                           key=lambda x: x[1], reverse=True))
        
        logger.info(f"\nEstimated serotype proportions:")
        for serotype, prop in estimated_proportions.items():
            in_mixture = "✓" if serotype in mixture_serotypes else "✗"
            logger.info(f"  {in_mixture} {serotype:20s}: {prop*100:6.2f}%")
        
        return estimated_proportions
    
    def evaluate_accuracy(self, true_serotypes, true_proportions, estimated_proportions):
        """
        Evaluate accuracy of serotype calling and quantification
        """
        logger.info("\n" + "=" * 80)
        logger.info("ACCURACY EVALUATION")
        logger.info("=" * 80)
        
        # Serotype detection
        true_set = set(true_serotypes)
        detected_set = set(estimated_proportions.keys())
        
        true_positives = len(true_set & detected_set)
        false_positives = len(detected_set - true_set)
        false_negatives = len(true_set - detected_set)
        
        sensitivity = true_positives / len(true_set) if len(true_set) > 0 else 0
        precision = true_positives / len(detected_set) if len(detected_set) > 0 else 0
        
        logger.info(f"\nSEROTYPE DETECTION:")
        logger.info(f"  True positives: {true_positives}")
        logger.info(f"  False positives: {false_positives}")
        logger.info(f"  False negatives: {false_negatives}")
        logger.info(f"  Sensitivity: {sensitivity*100:.1f}%")
        logger.info(f"  Precision: {precision*100:.1f}%")
        
        if false_positives > 0:
            logger.warning(f"\n  False positive serotypes detected:")
            for st in (detected_set - true_set):
                logger.warning(f"    {st}: {estimated_proportions[st]*100:.2f}%")
        
        if false_negatives > 0:
            logger.warning(f"\n  False negative serotypes (missed):")
            for st in (true_set - detected_set):
                logger.warning(f"    {st}")
        
        # Quantification accuracy (for correctly detected serotypes)
        logger.info(f"\nQUANTIFICATION ACCURACY:")
        
        errors = []
        for serotype, true_prop in zip(true_serotypes, true_proportions):
            if serotype in estimated_proportions:
                est_prop = estimated_proportions[serotype]
                absolute_error = abs(est_prop - true_prop)
                relative_error = absolute_error / true_prop if true_prop > 0 else 0
                
                errors.append(absolute_error)
                
                logger.info(f"  {serotype:20s}: True={true_prop*100:5.1f}%  Est={est_prop*100:5.1f}%  "
                          f"Error={absolute_error*100:5.1f}%")
        
        if errors:
            mae = np.mean(errors)
            logger.info(f"\n  Mean Absolute Error: {mae*100:.2f}%")
        
        return {
            'true_positives': true_positives,
            'false_positives': false_positives,
            'false_negatives': false_negatives,
            'sensitivity': sensitivity,
            'precision': precision,
            'mae': np.mean(errors) if errors else None
        }
    
    def run_simulation_suite(self, n_replicates=10):
        """
        Run a suite of mixture simulations
        """
        logger.info("\n" + "=" * 80)
        logger.info("RUNNING SIMULATION SUITE")
        logger.info("=" * 80)
        
        results = []
        
        # Define test mixtures
        test_mixtures = [
            # 2-serotype mixtures (equal proportions)
            (['19F', '19A'], [0.5, 0.5], 'equal_19F_19A'),
            (['6A', '6B'], [0.5, 0.5], 'equal_6A_6B'),
            (['23F', '23A'], [0.5, 0.5], 'equal_23F_23A'),
            
            # 2-serotype mixtures (unequal proportions)
            (['19F', '19A'], [0.8, 0.2], 'unequal_19F_19A_80_20'),
            (['19F', '19A'], [0.9, 0.1], 'unequal_19F_19A_90_10'),
            
            # 3-serotype mixtures
            (['19F', '19A', '6A'], [0.5, 0.3, 0.2], 'triple_19F_19A_6A'),
            (['14', '23F', '6B'], [0.4, 0.4, 0.2], 'triple_14_23F_6B'),
        ]
        
        for serotypes, proportions, name in test_mixtures:
            # Check if all serotypes exist
            if not all(st in self.reference_loci for st in serotypes):
                logger.warning(f"Skipping {name}: some serotypes not in database")
                continue
            
            logger.info(f"\n{'='*80}")
            logger.info(f"TEST: {name}")
            logger.info(f"{'='*80}")
            
            # Run replicates
            replicate_results = []
            
            for rep in range(n_replicates):
                logger.info(f"\nReplicate {rep+1}/{n_replicates}")
                
                # Generate mixture
                reads, read_map = self.simulate_mixture(serotypes, proportions, f"{name}_rep{rep}")
                
                # Map reads
                primer_counts = self.map_reads_to_amplicons(reads, read_map)
                
                # Deconvolve
                estimated = self.deconvolve_mixture(primer_counts, serotypes)
                
                # Evaluate
                accuracy = self.evaluate_accuracy(serotypes, proportions, estimated)
                
                replicate_results.append({
                    'mixture': name,
                    'replicate': rep,
                    'true_serotypes': ','.join(serotypes),
                    'true_proportions': ','.join([f"{p:.2f}" for p in proportions]),
                    **accuracy
                })
            
            results.extend(replicate_results)
        
        # Save results
        results_df = pd.DataFrame(results)
        output_file = self.output_dir / "results" / "simulation_results.csv"
        results_df.to_csv(output_file, index=False)
        
        logger.info(f"\n{'='*80}")
        logger.info(f"Saved results: {output_file}")
        
        return results_df
    
    def plot_results(self, results_df):
        """Plot simulation results"""
        logger.info("\nGenerating result plots...")
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Sensitivity by mixture type
        mixture_sensitivity = results_df.groupby('mixture')['sensitivity'].agg(['mean', 'std'])
        axes[0, 0].bar(range(len(mixture_sensitivity)), mixture_sensitivity['mean'], 
                      yerr=mixture_sensitivity['std'], capsize=5, alpha=0.7)
        axes[0, 0].set_xticks(range(len(mixture_sensitivity)))
        axes[0, 0].set_xticklabels(mixture_sensitivity.index, rotation=45, ha='right', fontsize=8)
        axes[0, 0].set_ylabel('Sensitivity', fontweight='bold')
        axes[0, 0].set_title('Serotype Detection Sensitivity', fontweight='bold')
        axes[0, 0].set_ylim(0, 1.1)
        axes[0, 0].axhline(y=1.0, color='r', linestyle='--', alpha=0.5)
        axes[0, 0].grid(True, alpha=0.3, axis='y')
        
        # Precision by mixture type
        mixture_precision = results_df.groupby('mixture')['precision'].agg(['mean', 'std'])
        axes[0, 1].bar(range(len(mixture_precision)), mixture_precision['mean'],
                      yerr=mixture_precision['std'], capsize=5, alpha=0.7, color='green')
        axes[0, 1].set_xticks(range(len(mixture_precision)))
        axes[0, 1].set_xticklabels(mixture_precision.index, rotation=45, ha='right', fontsize=8)
        axes[0, 1].set_ylabel('Precision', fontweight='bold')
        axes[0, 1].set_title('Serotype Detection Precision', fontweight='bold')
        axes[0, 1].set_ylim(0, 1.1)
        axes[0, 1].axhline(y=1.0, color='r', linestyle='--', alpha=0.5)
        axes[0, 1].grid(True, alpha=0.3, axis='y')
        
        # MAE by mixture type
        mixture_mae = results_df.groupby('mixture')['mae'].agg(['mean', 'std'])
        axes[1, 0].bar(range(len(mixture_mae)), mixture_mae['mean']*100,
                      yerr=mixture_mae['std']*100, capsize=5, alpha=0.7, color='orange')
        axes[1, 0].set_xticks(range(len(mixture_mae)))
        axes[1, 0].set_xticklabels(mixture_mae.index, rotation=45, ha='right', fontsize=8)
        axes[1, 0].set_ylabel('Mean Absolute Error (%)', fontweight='bold')
        axes[1, 0].set_title('Quantification Accuracy', fontweight='bold')
        axes[1, 0].grid(True, alpha=0.3, axis='y')
        
        # False positive rate
        fp_rate = results_df.groupby('mixture')['false_positives'].mean()
        axes[1, 1].bar(range(len(fp_rate)), fp_rate, alpha=0.7, color='red')
        axes[1, 1].set_xticks(range(len(fp_rate)))
        axes[1, 1].set_xticklabels(fp_rate.index, rotation=45, ha='right', fontsize=8)
        axes[1, 1].set_ylabel('False Positives', fontweight='bold')
        axes[1, 1].set_title('False Positive Rate', fontweight='bold')
        axes[1, 1].grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        output_file = self.output_dir / "plots" / "simulation_results.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved: {output_file}")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Simulate mixed samples and validate primer performance'
    )
    
    parser.add_argument('--database', required=True)
    parser.add_argument('--primers', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--coverage', type=int, default=100,
                       help='Sequencing coverage depth (default: 100)')
    parser.add_argument('--replicates', type=int, default=10,
                       help='Number of replicates per mixture (default: 10)')
    parser.add_argument('--max-mismatches', type=int, default=2)
    parser.add_argument('--debug', action='store_true')
    
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    simulator = MixtureSimulator(
        database_dir=args.database,
        primer_file=args.primers,
        output_dir=args.output,
        max_mismatches=args.max_mismatches,
        coverage_depth=args.coverage
    )
    
    try:
        logger.info("=" * 80)
        logger.info("IN SILICO MIXTURE VALIDATION")
        logger.info("=" * 80)
        
        simulator.load_reference_loci()
        simulator.load_primers()
        simulator.map_primers()
        
        results_df = simulator.run_simulation_suite(n_replicates=args.replicates)
        simulator.plot_results(results_df)
        
        logger.info("\n" + "=" * 80)
        logger.info("VALIDATION COMPLETE!")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Failed: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()