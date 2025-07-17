#!/bin/bash

# SLURM directives for job submission
# -n 4: Request 4 tasks (cores)
# -p Cnode_all: Specify the partition (queue) to use
#SBATCH -n 4
#SBATCH -p Cnode_all

# Script Purpose:
# This script automates peak calling using MACS3 on a set of BAM files in a SLURM cluster.
# It processes each BAM file in the input directory, using a common control BAM file.
# Outputs are saved to the specified output directory, with logging and checks for existing files.

# Set environment variables and parameters
CONTROL="/public/home/huo/DAP/bam/Control/ck.dedup.sorted.filter.bam"  # Path to the control BAM file
INPUT_DIR="/public/home/huo/DAP/bam/"  # Directory containing input BAM files
OUTPUT_DIR="/public/home/huo/DAP/macs3_results/"  # Directory for output files (peaks and logs)
GENOME_SIZE="1804145610"  # Effective genome size for MACS3 (e.g., for maize or custom genome)

# Ensure the output directory exists; create it if not
mkdir -p "$OUTPUT_DIR"

# Check if the control file exists; exit if not
if [ ! -f "$CONTROL" ]; then
    echo "Error: Control file not found: $CONTROL"
    exit 1
fi

# Check if the input directory exists; exit if not
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Create a subdirectory for logs within the output directory
mkdir -p "$OUTPUT_DIR/logs"

# Process each BAM file in the input directory
# Pattern: Matches files ending with "_processed.filter.sorted.dedup.bam"
for bam_file in "$INPUT_DIR"/*_processed.filter.sorted.dedup.bam; do
    if [ -f "$bam_file" ]; then
        # Extract the base name by removing the suffix from the file name
        base_name=$(basename "$bam_file" _processed.filter.sorted.dedup.bam)
        
        # Define the expected output file (narrowPeak format)
        output_file="$OUTPUT_DIR/${base_name}_peaks.narrowPeak"
        
        # Check if the output file already exists; skip processing if it does
        if [ -f "$output_file" ]; then
            echo "Output file already exists for: $base_name. Skipping processing."
            continue
        fi
        
        # Print status message
        echo "Processing: $base_name"
        
        # Run MACS3 callpeak command
        # -t: Treatment BAM file
        # -c: Control BAM file
        # -n: Output prefix
        # -f BAMPE: Input format is BAM paired-end
        # --call-summits: Call summits within peaks
        # --seed: Random seed for reproducibility
        # --outdir: Output directory
        # -g: Genome size
        # -q: q-value cutoff (0.05)
        # --keep-dup all: Keep all duplicates
        # Redirect stdout and stderr to a log file
        macs3 callpeak \
            -t "$bam_file" \
            -c "$CONTROL" \
            -n "${base_name}" \
            -f BAMPE \
            --call-summits \
            --seed 666 \
            --outdir "$OUTPUT_DIR" \
            -g "$GENOME_SIZE" \
            -q 0.05 \
            --keep-dup all \
            > "$OUTPUT_DIR/logs/${base_name}_macs3.log" 2>&1
        
        # Check if the MACS3 command succeeded (exit code 0)
        if [ $? -eq 0 ]; then
            echo "Completed processing: $base_name"
        else
            echo "Error processing: $base_name. Check log for details."
        fi
    fi
done

# Final status message after all processing
echo "All peak calling completed!"
