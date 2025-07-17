#!/bin/bash

# Optional SLURM directives for job submission (comment out if not needed)
# -n 4: Request 4 tasks (cores) for parallel processing
# -p Cnode: Specify the partition (queue) to use
# SBATCH -n 4
# SBATCH -p Cnode

# Script Purpose:
# This script runs MSPC (Multi-Sample Peak Caller) on pairs of BED files (e.g., replicates named *-1.bed and *-2.bed)
# from an input directory. It generates consensus peaks in the output directory, skips pairs if outputs exist,
# and handles errors. Assumes MSPC is installed and the executable path is provided.

# Set locale to ensure consistent behavior (e.g., for sorting and numeric formats)
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# Define input and output directories
INPUT_DIR="/public/home/huo/DAP/macs3"  # Directory containing input BED files
OUTPUT_DIR="/public/home/huo/DAP/MSPC"  # Directory for output consensus peaks

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define the absolute path to the MSPC executable
MSPC_EXEC="/public/home/huo/DAP/MSPC"

# Check if the MSPC executable exists; exit if not
if [ ! -x "$MSPC_EXEC" ]; then
    echo "Error: MSPC executable not found or not executable: $MSPC_EXEC"
    exit 1
fi

# Check if input directory exists; exit if not
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Loop over all BED files ending with -1.bed (assuming pairs like *-1.bed and *-2.bed)
for file1 in "$INPUT_DIR"/*-1.bed; do
    # Extract base name by removing the -1.bed suffix
    base_name=$(basename "$file1" -1.bed)
    
    # Define the paired file (replace -1 with -2)
    file2="$INPUT_DIR/${base_name}-2.bed"
    
    # Define the output file path (ConsensusPeaks.bed in a subdirectory)
    output_file="$OUTPUT_DIR/$base_name/ConsensusPeaks.bed"

    # Check if the paired file exists; skip if not
    if [ ! -f "$file2" ]; then
        echo "Warning: Paired file not found: $file2"
        continue
    fi

    # Check if the output file already exists; skip if it does
    if [ -f "$output_file" ]; then
        echo "Skipping: Output file already exists: $output_file"
        continue
    fi

    # Run MSPC command
    # -i: Input BED files (multiple can be provided)
    # -r technical: Replication type (technical replicates)
    # -w 1e-2: Weak threshold for peak calling
    # -s 1e-4: Stringency threshold
    # -a 0.05: Alpha (significance level)
    # -c 2: Minimum number of replicates required
    # -m lowest: Method for combining p-values (lowest)
    # -d 64: Tau value (distance parameter)
    # -o: Output directory
    # --excludeHeader true: Exclude headers from output
    if "$MSPC_EXEC" -i "$file1" "$file2" -r technical -w 1e-2 -s 1e-4 -a 0.05 -c 2 -m lowest -d 64 -o "$OUTPUT_DIR/$base_name" --excludeHeader true; then
        echo "Processing completed: $output_file"
    else
        echo "Error: Failed to process files $file1 and $file2"
    fi
done

# Final status message
echo "All BED file pairs processed!"
