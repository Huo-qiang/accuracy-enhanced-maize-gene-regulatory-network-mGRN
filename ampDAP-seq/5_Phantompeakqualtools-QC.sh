#!/bin/bash

# SLURM directives for job submission
# -n 32: Request 32 tasks (cores)
# -p Cnode_all: Specify the partition (queue) to use
# -N 1: Request 1 node
#SBATCH -n 32
#SBATCH -p Cnode_all
#SBATCH -N 1

# Script Purpose:
# This script runs quality control (QC) on BAM files using the run_spp.R script from phantompeakqualtools
# in a SLURM cluster environment. It processes each BAM file concurrently (limited by max_jobs),
# generates PDF plots (-savp) and text outputs (-out), and skips files if outputs already exist.
# Assumes phantompeakqualtools is installed in the specified conda environment.

# Activate the conda environment using mamba
mamba activate phantompeakqualtools_env

# Set the number of threads for OpenBLAS to avoid overloading (e.g., for R computations)
export OPENBLAS_NUM_THREADS=10

# Maximum number of concurrent background jobs
max_jobs=10
job_count=0  # Counter for current running jobs

# Loop over all BAM files in the specified directory (pattern: *_processed.filter.sorted.dedup.bam)
for sample in /public/home/huo/DAP/bam/*_processed.filter.sorted.dedup.bam; do
    # Extract base name by removing the suffix
    base=$(basename "$sample" "_processed.filter.sorted.dedup.bam")

    # Define output file paths
    pdf_output="/public/home/huo/bam_qc/${base}.pdf"
    txt_output="/public/home/huo/bam_qc/${base}.txt"

    # Check if both output files already exist; skip if they do
    if [[ -f $pdf_output && -f $txt_output ]]; then
        echo "results ${base} exits，skip..."
        continue
    fi

    # If outputs do not exist, run Rscript in the background
    echo "proccessed ${base}..."
    Rscript run_spp.R \
        -c="$sample" \  # Input BAM file
        -p=10 \  # Number of processors (threads) for run_spp.R
        -savp="$pdf_output" \  # Save plot to PDF
        -out="$txt_output" &  # Save output to text file

    # Increment the job counter
    ((job_count++))

    # If the maximum concurrent jobs is reached, wait for them to finish
    if (( job_count >= max_jobs )); then
        # Wait for all background jobs to complete
        wait
        # Reset the job counter after waiting
        job_count=0
    fi
done

# Wait for any remaining background jobs to complete
wait

# Final status message (optional: add to indicate script completion)
echo "All BAM QC processing completed!"
