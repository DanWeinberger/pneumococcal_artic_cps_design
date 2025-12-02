#!/bin/bash
#SBATCH --job-name=optimize_primers
#SBATCH --output=optimize_primers_%j.log
#SBATCH --error=optimize_primers_%j.err
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=day

# Primer Set Optimization Script for SLURM
# This script optimizes a large primer set for serotype discrimination

echo "========================================"
echo "Primer Set Optimization"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# ============================================================================
# CONFIGURATION - EDIT THESE PATHS
# ============================================================================

# Path to your database directory (containing reference.fasta files)
DATABASE_DIR="/gpfs/gibbs/project/weinberger_daniel/dmw63/artic_serotype/seroba/database"

# Path to your input primers file (from pangraph_primer_design.py)
PRIMERS_FILE="/gpfs/gibbs/project/weinberger_daniel/dmw63/artic_serotype/seroba/pangraph_output/primers/pangraph_primers.csv"

# Output directory for optimized primers
OUTPUT_DIR="/gpfs/gibbs/project/weinberger_daniel/dmw63/artic_serotype/seroba/optimized_primers"

# Path to the optimization script
SCRIPT_PATH="/gpfs/gibbs/project/weinberger_daniel/dmw63/artic_serotype/seroba/optimize_primer_set.py"

# ============================================================================
# OPTIMIZATION PARAMETERS
# ============================================================================

# Target amplicons per serotype
MIN_AMPLICONS=5
MAX_AMPLICONS=10

# Maximum mismatches allowed for primer binding
MAX_MISMATCHES=2

# Output format (csv or tsv)
FORMAT="csv"

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

# Load required modules (adjust for your cluster)
# module load Python/3.9
# module load Anaconda3

# If using conda environment, activate it
# source activate primer_design

# Or if using your base conda
# eval "$(conda shell.bash hook)"
# conda activate base

echo "Python version:"
python --version
echo ""

echo "Checking required Python packages..."
python -c "import pandas; import numpy; import Bio; print('✓ All required packages available')"
echo ""

# ============================================================================
# RUN OPTIMIZATION
# ============================================================================

echo "========================================"
echo "Running Primer Optimization"
echo "========================================"
echo "Database: $DATABASE_DIR"
echo "Input primers: $PRIMERS_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Target: $MIN_AMPLICONS-$MAX_AMPLICONS amplicons per serotype"
echo ""

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run the optimization
python "$SCRIPT_PATH" \
    --database "$DATABASE_DIR" \
    --primers "$PRIMERS_FILE" \
    --output "$OUTPUT_DIR" \
    --min-amplicons "$MIN_AMPLICONS" \
    --max-amplicons "$MAX_AMPLICONS" \
    --max-mismatches "$MAX_MISMATCHES" \
    --format "$FORMAT"

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "Optimization completed successfully!"
    echo "========================================"
    echo ""
    echo "Results saved to: $OUTPUT_DIR"
    echo ""
    echo "Output files:"
    ls -lh "$OUTPUT_DIR/optimized_primers/"
    echo ""
    echo "Report:"
    cat "$OUTPUT_DIR/reports/optimization_report.txt"
else
    echo ""
    echo "========================================"
    echo "ERROR: Optimization failed!"
    echo "========================================"
    exit 1
fi

echo ""
echo "End time: $(date)"
echo "========================================"