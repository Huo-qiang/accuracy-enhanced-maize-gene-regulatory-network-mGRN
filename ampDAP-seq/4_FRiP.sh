#!/bin/bash

# Optional SLURM directives for job submission (comment out if not needed)
# -p Cnode: Specify the partition (queue) to use
# -n 4: Request 4 tasks (cores) for parallel processing
# SBATCH -p Cnode
# SBATCH -n 4

# Script Purpose:
# This script calculates the Fraction of Reads in Peaks (FRiP) score for each BAM file in a directory.
# For each BAM file:
#   1. Counts total reads using Samtools.
#   2. Counts reads overlapping with peaks (from macs3 results) using Bedtools.
#   3. Computes FRiP as (reads_in_peaks / total_reads).
#   4. Outputs the sample name, FRiP score, and label to a file in the FRiP directory.
# Assumes peaks files are in BED format and named consistently with BAM files.

# Set directory paths
dirdedupbam="/public/home/huo/DAP/bam"  # Directory containing deduplicated sorted BAM files
dirpeak="/public/home/huo/DAP/macs3_results"  # Directory containing macs3 peak files 
dirFRIP="/public/home/huo/DAP/FRiP"  # Directory for output FRiP score files

# Create output directory if it doesn't exist
mkdir -p "$dirFRIP"

# Check if required directories exist; exit if not
for dir in "$dirdedupbam" "$dirgem" "$dirFRIP"; do
    if [ ! -d "$dir" ]; then
        echo "Error: Directory not found: $dir"
        exit 1
    fi
done

# Loop over all BAM files (pattern: *.dedup.sorted.bam)
for sample in "$dirdedupbam"/*.dedup.sorted.bam; do
    # Extract base name by removing the suffix
    base=$(basename "$sample" ".dedup.sorted.bam")
    echo "Processing sample: ${base}"

    # Check if the BAM file exists (safety check)
    if [ ! -f "${dirdedupbam}/${base}.dedup.sorted.bam" ]; then
        echo "Error: BAM file not found for sample: ${base}"
        continue
    fi

    # Check if the corresponding peaks file exists; skip if not
    if [ ! -f "${dirgem}/${base}.GEM_rmgreylist.bed" ]; then
        echo "Error: Peaks file not found for sample: ${base}"
        continue
    fi

    # Calculate total reads in the BAM file
    total_reads=$(samtools view -c "${dirdedupbam}/${base}.dedup.sorted.bam")
    if [ $? -ne 0 ]; then
        echo "Error: Failed to count total reads for sample: ${base}"
        continue
    fi

    # Calculate reads in peaks:
    # - Sort the peaks BED file
    # - Merge overlapping peaks
    # - Intersect with BAM (unique reads only, output as UBAM)
    # - Count the resulting reads
    reads_in_peaks=$(bedtools sort -i "${dirgem}/${base}.GEM_rmgreylist.bed" \
        | bedtools merge -i stdin \
        | bedtools intersect -u -nonamecheck \
            -a "${dirdedupbam}/${base}.dedup.sorted.bam" \
            -b stdin -ubam \
        | samtools view -c)
    if [ $? -ne 0 ]; then
        echo "Error: Failed to count reads in peaks for sample: ${base}"
        continue
    fi

    # Calculate FRiP score using awk (handle potential division by zero implicitly)
    frip_score=$(awk "BEGIN {print (${reads_in_peaks} / ${total_reads})}")

    # Output to file: sample name, FRiP score, and label "FRiP"
    echo "${base}" > "${dirFRIP}/${base}.FRiP"
    echo "${frip_score}" >> "${dirFRIP}/${base}.FRiP"
    echo "FRiP" >> "${dirFRIP}/${base}.FRiP"

    echo "Completed processing sample: ${base} (FRiP: ${frip_score})"
done

# Final status message
echo "All samples processed successfully!"
