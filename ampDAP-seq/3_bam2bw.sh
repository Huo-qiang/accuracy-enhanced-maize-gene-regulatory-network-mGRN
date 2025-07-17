#!/bin/bash

# SLURM directives for job submission
# -p Cnode: Specify the partition (queue) to use
# -n 64: Request 64 tasks (cores)
#SBATCH -p Cnode
#SBATCH -n 64

# Script Purpose:
# This script converts BAM files to BigWig (.bw) files using bamCoverage from deepTools in a SLURM cluster.
# It processes each BAM file in the input directory, normalizes using RPKM, and outputs to the specified directory.
# Includes concurrency control (limited to max_jobs), skips existing outputs, and cleans up with wait commands.

# Set paths
dirbam="/public/home/huo/bam"  # Directory containing input BAM files
dirbw="/public/home/huo/bw"  # Directory for output BigWig files

# Create output directory if it doesn't exist
mkdir -p "$dirbw"

# Maximum number of concurrent background jobs
max_jobs=5
job_count=0  # Counter for current running jobs

# Loop over all BAM files (pattern: *.filter.sorted.dedup.bam)
for sample in "$dirbam"/*.filter.sorted.dedup.bam; do
    # Extract base name by removing the suffix
    base=$(basename "$sample" ".filter.sorted.dedup.bam")
    output_file="$dirbw/${base}.bw"

    # Check if the output file already exists; skip if it does
    if [[ -f "$output_file" ]]; then
        echo "Output file $output_file already exists. Skipping..."
        continue
    fi

    # Run bamCoverage to generate BigWig file
    # -b: Input BAM file
    # -o: Output BigWig file
    # --binSize 10: Bin size for coverage calculation
    # --normalizeUsing RPKM: Normalize using Reads Per Kilobase Million
    # -p 64: Number of processors (threads) to use
    bamCoverage -b "$sample" -o "$output_file" --binSize 10 --normalizeUsing RPKM -p 64 &

    # Increment the job counter
    ((job_count++))

    # If the maximum concurrent jobs is reached, wait for them to finish
    if [[ $job_count -ge $max_jobs ]]; then
        wait  # Wait for all background jobs to complete
        job_count=0  # Reset the job counter
    fi
done

# Wait for any remaining background jobs to complete
wait

# Final status message
echo "All processing completed!"
