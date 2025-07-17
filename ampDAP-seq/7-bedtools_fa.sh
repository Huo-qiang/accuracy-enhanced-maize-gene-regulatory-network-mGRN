#!/bin/bash

# Optional SLURM directives for job submission (comment out if not needed)
# -n 64: Request 64 tasks (cores) for high concurrency
# -p Cnode_all: Specify the partition (queue) to use
# -N 1: Request 1 node
# SBATCH -n 64
# SBATCH -p Cnode_all
# SBATCH -N 1

# Script Purpose:
# This script processes BED files from peak calling (e.g., MACS3 outputs) to extract FASTA sequences.
# For each BED file:
#   1. Selects the top 600 peaks.
#   2. Calculates the center of each peak and extends it by 100bp on each side (total 200bp).
#   3. Uses bedtools getfasta to extract sequences from a reference genome.
# It runs processes concurrently (limited by max_jobs) for efficiency, using background jobs and temporary files.

# Define input and output directories
dirmeme="/mnt/disk16t/DAP/new_DAP/macs3/161_newck_q0.05/resized_peaks/sort"  # Directory containing input BED files (*_filtered.bed)
dirfa="/mnt/disk16t/DAP/new_DAP/fa_161_top600_200bp_before"  # Directory for output FASTA files

# Set maximum concurrent jobs
max_jobs=100
job_count=0  # Counter for current running jobs

# Define reference genome FASTA file
ref_fasta="/mnt/disk16t/DAP/new_DAP/Zea_mays.B73_RefGen_v4.49.dna.chr.fa"

# Ensure output directory exists
mkdir -p "$dirfa"

# Check if required tools and files exist; exit if not
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is not installed or not in PATH"
    exit 1
fi
if ! command -v awk &> /dev/null; then
    echo "Error: awk is not installed or not in PATH"
    exit 1
fi
if [ ! -f "$ref_fasta" ]; then
    echo "Error: Reference FASTA file not found: $ref_fasta"
    exit 1
fi
if [ ! -d "$dirmeme" ]; then
    echo "Error: Input directory not found: $dirmeme"
    exit 1
fi

# Function to wait for available job slots (reduces to below max_jobs)
wait_for_jobs() {
    while [ $job_count -ge $max_jobs ]; do
        wait -n  # Wait for one background job to finish
        ((job_count--))  # Decrement the job counter
    done
}

# Loop over all input BED files (pattern: *_filtered.bed)
for sample in "$dirmeme"/*_filtered.bed; do
    # Check if the file exists (skip if no files match the pattern)
    if [ ! -f "$sample" ]; then
        echo "Warning: No *_filtered.bed files found in $dirmeme"
        continue
    fi

    # Extract base name by removing the _filtered.bed suffix
    base=$(basename "$sample" "_filtered.bed")
    
    # Wait for an available job slot before starting a new one
    wait_for_jobs
    
    # Run processing in a background subshell
    (
        # Create a temporary BED file for processed peaks
        temp_bed=$(mktemp)
        
        # Process the top 600 peaks:
        # - Select first 600 lines (NR <= 600)
        # - Calculate center as integer average of start ($2) and end ($3)
        # - Extend 100bp on each side (start = center - 100, end = center + 100)
        # - Ensure start >= 0
        # - Output: chr \t start \t end
        awk 'NR <= 600 {
            center = int(($2 + $3)/2);
            start = center - 100;
            end = center + 100;
            if (start < 0) start = 0;
            print $1 "\t" start "\t" end;
        }' "$sample" > "$temp_bed"
        
        # Check if awk succeeded
        if [ $? -ne 0 ]; then
            echo "Error: awk processing failed for $base"
            rm "$temp_bed"
            exit 1
        fi
        
        # Extract FASTA sequences using bedtools getfasta
        # -fi: Reference FASTA file
        # -bed: Input BED file (processed temp_bed)
        # -fo: Output FASTA file
        bedtools getfasta -fi "$ref_fasta" \
                          -bed "$temp_bed" \
                          -fo "${dirfa}/${base}.fa"
        
        # Check if bedtools succeeded
        if [ $? -ne 0 ]; then
            echo "Error: bedtools getfasta failed for $base"
            rm "$temp_bed"
            exit 1
        fi
        
        # Clean up temporary file
        rm "$temp_bed"
        
        echo "Completed processing: $base"
    ) &
    
    # Increment the job counter
    ((job_count++))
    
    echo "Started processing: $base (Current jobs: $job_count)"
done

# Wait for all remaining background jobs to complete
wait

# Final status message
echo "All files have been processed"
