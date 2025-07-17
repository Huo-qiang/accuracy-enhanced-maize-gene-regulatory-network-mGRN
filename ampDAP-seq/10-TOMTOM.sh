#!/bin/bash

# Optional SLURM directives for job submission (comment out if not needed)
# -p amd_256: Specify the partition (queue) to use
# -N 1: Request 1 node
# -n 1: Request 1 task (core) – suitable for single-threaded jobs
# SBATCH -p amd_256
# SBATCH -N 1
# SBATCH -n 1

# Script Purpose:
# This script runs TomTom (from MEME suite) to compare motifs in input .meme files (e.g., from MEME-ChIP)
# against two reference motif databases. It identifies similarities between query motifs and known motifs,
# processes each file sequentially, and outputs results to subdirectories. Skips files if outputs already exist.

# Define paths
dirmeme="/public/home/huo/DAP/memechip/meme"  # Input directory containing .meme files
dirout="/public/home/huo/DAP/tomtom"  # Output directory for TomTom results
reference1="JASPAR2024_CORE_non-redundant_pfms.meme"   # Arabidopsis dap-seq motif database
reference2="ArabidopsisDAPv1.meme"  # jaspar2024 motif database

# Define the TomTom executable path
tomtom_executable="/public/home/huo/meme-5.5.1/src/tomtom"

# Ensure output directory exists
mkdir -p "$dirout"

# Check if required directories, references, and tools exist; exit if not
if [ ! -d "$dirmeme" ]; then
    echo "Error: Input directory not found: $dirmeme"
    exit 1
fi
if [ ! -f "$reference1" ]; then
    echo "Error: Reference database 1 not found: $reference1"
    exit 1
fi
if [ ! -f "$reference2" ]; then
    echo "Error: Reference database 2 not found: $reference2"
    exit 1
fi
if [ ! -x "$tomtom_executable" ]; then
    echo "Error: TomTom executable not found or not executable: $tomtom_executable"
    exit 1
fi

# Loop over all .meme files (pattern: *.meme)
for sample in "$dirmeme"/*.meme; do
    # Check if the file exists (skip if no files match the pattern)
    if [ ! -f "$sample" ]; then
        echo "Warning: No *.meme files found in $dirmeme"
        continue
    fi

    # Extract base name by removing the .meme suffix
    base=$(basename "$sample" ".meme")

    # Define output subdirectory
    output_subdir="${dirout}/${base}"

    # Check if the output subdirectory already exists; skip if it does
    if [ -d "$output_subdir" ]; then
        echo "Skipping: Output directory already exists for $base: $output_subdir"
        continue
    fi

    # Run TomTom
    # -oc: Output directory (creates a subdirectory)
    # -min-overlap 5: Minimum overlap between motifs
    # -dist pearson: Use Pearson correlation as the distance metric
    # -evalue: Report E-values
    # -thresh 1: E-value threshold for significance
    # -no-ssc: Do not use small-sample correction
    # Input: Query .meme file, followed by reference databases
    if "$tomtom_executable" -oc "$output_subdir" -min-overlap 5 -dist pearson -evalue -thresh 1 -no-ssc "$sample" "$reference1" "$reference2"; then
        echo "Processed $sample -> $output_subdir"
    else
        echo "Error: TomTom failed for $base"
    fi
done

# Final status message
echo "All .meme files have been processed!"
