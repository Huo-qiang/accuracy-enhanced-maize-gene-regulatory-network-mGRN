#!/bin/bash

# Optional SLURM directives for job submission (comment out if not needed)
# -p amd_256: Specify the partition (queue) to use
# -N 1: Request 1 node
# -n 1: Request 1 task (core) – adjust if higher parallelism is needed
# SBATCH -p amd_256
# SBATCH -N 1
# SBATCH -n 1

# Script Purpose:
# This script runs AME (Analysis of Motif Enrichment) from the MEME suite on FASTA files for motif enrichment analysis.
# It uses two motif databases, processes files concurrently (limited by max_jobs), and outputs results to subdirectories.
# Skips files if outputs already exist, shuffles sequences for control, and applies specified scoring/thresholds.

# Define paths
input_dir="/public/home/huo/DAP/fa"  # Input directory containing FASTA files (*.fa)
output_dir="/public/home/huo/DAP/AME"  # Output directory for AME results
db1="/mnt/disk16t/DAP/new_DAP/ArabidopsisDAPv1.meme"  # Arabidopsis dap-seq motif database
db2="/mnt/disk16t/DAP/new_DAP/JASPAR2024_CORE_non-redundant_pfms.meme"  # jaspar2024 motif database

# Define the AME executable path
ame_exec="/public/home/meme-5.5.1/src/ame"

# Maximum number of concurrent background jobs
max_jobs=10
job_count=0  # Counter for current running jobs

# Ensure output directory exists
mkdir -p "$output_dir"

# Check if required directories, databases, and tools exist; exit if not
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory not found: $input_dir"
    exit 1
fi
if [ ! -f "$db1" ]; then
    echo "Error: Database 1 not found: $db1"
    exit 1
fi
if [ ! -f "$db2" ]; then
    echo "Error: Database 2 not found: $db2"
    exit 1
fi
if [ ! -x "$ame_exec" ]; then
    echo "Error: AME executable not found or not executable: $ame_exec"
    exit 1
fi

# Loop over all FASTA files (pattern: *.fa)
for sample in "$input_dir"/*.fa; do
    # Check if the input file exists; skip if not
    if [ ! -f "$sample" ]; then
        echo "Warning: Input file not found: $sample"
        continue
    fi

    # Extract base name by removing the .fa suffix
    base=$(basename "$sample" ".fa")

    # Define output subdirectory
    output_subdir="${output_dir}/${base}"

    # Check if the output subdirectory already exists; skip if it does
    if [ -d "$output_subdir" ]; then
        echo "Skipping: Output directory already exists for $base: $output_subdir"
        continue
    fi

    # Run AME in the background
    # --verbose 2: Set verbosity level
    # --oc: Output directory (creates a subdirectory)
    # --scoring avg: Use average scoring method
    # --method fisher: Use Fisher's exact test
    # --hit-lo-fraction 0.25: Fraction of low-scoring sites to consider hits
    # --evalue-report-threshold 0.05: E-value threshold for reporting
    # --control --shuffle--: Use shuffled sequences as control
    # --kmer 2: K-mer length for shuffling
    # Input: FASTA file, followed by motif databases
    "$ame_exec" --verbose 2 --oc "$output_subdir" --scoring avg \
        --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 0.05 \
        --control --shuffle-- --kmer 2 "$sample" "$db1" "$db2" &

    # Increment the job counter
    job_count=$((job_count + 1))

    # If the maximum concurrent jobs is reached, wait for them to finish
    if [ "$job_count" -ge "$max_jobs" ]; then
        wait  # Wait for all background jobs to complete
        job_count=0  # Reset the job counter
    fi

    echo "Started processing: $base"
done

# Wait for any remaining background jobs to complete
wait

# Final status message
echo "All processing completed!"
