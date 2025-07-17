#!/bin/bash

# SLURM directives for job submission
# -p amd_256: Specify the partition (queue) to use
# -N 1: Request 1 node
# -n 1: Request 1 task (core) – suitable for single-threaded or low-parallel jobs
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 1

# Script Purpose:
# This script runs MEME-ChIP for motif discovery on FASTA files (e.g., sequences extracted from peaks).
# It processes each FASTA file in the input directory, discovers motifs with specified parameters
# (e.g., width 6-12, 30 processors for MEME), skips SpaMo and FIMO, and outputs to subdirectories.
# Skips processing if the output directory already exists.

# Define paths
dirout="/public/home/huo/DAP/fa/memechip"  # Output directory for MEME-ChIP results
dirfa="/public/home/huo/DAP/fa"  # Input directory containing FASTA files (*.fa)
# dirdb="/mnt/disk16t/DAP/GEM/JASPAR2024_CORE_non-redundant_pfms.meme"  # (Commented out in original; uncomment if needed for custom database)

# Define the meme-chip executable path
meme_chip_exec="/public/home/meme-5.5.1/scripts/meme-chip"

# Ensure output directory exists
mkdir -p "$dirout"

# Check if required directories and tools exist; exit if not
if [ ! -d "$dirfa" ]; then
    echo "Error: Input directory not found: $dirfa"
    exit 1
fi
if [ ! -x "$meme_chip_exec" ]; then
    echo "Error: meme-chip executable not found or not executable: $meme_chip_exec"
    exit 1
fi

# Loop over all FASTA files (pattern: *.fa)
for sample in "$dirfa"/*.fa; do
    # Check if the file exists (skip if no files match the pattern)
    if [ ! -f "$sample" ]; then
        echo "Warning: No *.fa files found in $dirfa"
        continue
    fi

    # Extract base name by removing the .fa suffix
    base=$(basename "$sample" ".fa")

    # Define output subdirectory
    output_subdir="${dirout}/${base}"

    # Check if the output subdirectory already exists; skip if it does
    if [ -d "$output_subdir" ]; then
        echo "Skipping: Output directory already exists for $base: $output_subdir"
        continue
    fi

    # Run meme-chip
    # -minw 6 -maxw 12: Motif width range (6 to 12)
    # -meme-p 30: Use 30 processors for MEME
    # -meme-nmotifs 0: Discover unlimited motifs (0 means no limit)
    # -meme-searchsize 0: No limit on search size
    # -spamo-skip: Skip SpaMo analysis
    # -fimo-skip: Skip FIMO analysis
    # -centrimo-local: Run CentriMo in local mode
    # -oc: Output directory (creates a subdirectory)
    if "$meme_chip_exec" -minw 6 -maxw 12 -meme-p 30 -meme-nmotifs 0 -meme-searchsize 0 -spamo-skip -fimo-skip -centrimo-local -oc "$output_subdir" "$sample"; then
        echo "Processed $sample -> $output_subdir"
    else
        echo "Error: meme-chip failed for $base"
    fi
done

# Final status message
echo "All FASTA files have been processed!"
