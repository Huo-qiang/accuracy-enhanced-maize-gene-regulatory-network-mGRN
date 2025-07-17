#!/bin/bash

# SLURM directives for job submission
# -n 64: Request 64 tasks (cores)
# -p Cnode_all: Specify the partition (queue) to use
# -N 1: Request 1 node
#SBATCH -n 64
#SBATCH -p Cnode_all
#SBATCH -N 1

# Script Purpose:
# This script processes paired-end FASTQ files to generate filtered, sorted, and deduplicated BAM files in a SLURM cluster.
# It uses Bowtie2 for alignment, Samtools for conversion/filtering/sorting/stats, and Sambamba for duplicate marking.
# Includes error handling, concurrency control (limited to max_jobs), and cleanup of intermediate files.
# If all steps succeed, it optionally deletes original FASTQ files (currently commented out).

# Exit on error (-e) and treat unset variables as errors (-u)
set -e
set -u

# Define main directories and paths
dirfastq="/public/home/huo/DAP/fq"  # Directory containing input FASTQ files
dirsam="/public/home/huo/DAP/bam"  # Directory for output SAM/BAM files
ref_index="/public/home/huo/b73v4_bowtie2/maizeb73v4_4.49"  # Path to Bowtie2 reference genome index (prefix)

# Check if necessary directories exist; exit if not
for dir in "$dirfastq" "$dirsam"; do
    if [ ! -d "$dir" ]; then
        echo "Error: Directory not found: $dir"
        exit 1
    fi
done

# Check if the reference genome index exists (check for .1.bt2 file as indicator)
if [ ! -e "${ref_index}.1.bt2" ]; then
    echo "Error: Bowtie2 index not found at $ref_index"
    exit 1
fi

# Set thread counts and other parameters
BOWTIE_THREADS=64  # Threads for Bowtie2
SAMTOOLS_THREADS=20  # Threads for Samtools
SAMBAMBA_THREADS=10  # Threads for Sambamba
QUALITY_THRESHOLD=30  # Minimum mapping quality for filtering
FILTER_FLAG=4  # SAM flag to filter out unmapped reads
max_jobs=1  # Maximum number of concurrent background jobs
job_count=0  # Counter for current running jobs

# Loop over each forward FASTQ file (pattern: *_1.fq.gz)
for sample in "${dirfastq}"/*_1.fq.gz; do
    # Extract base name by removing suffix
    base=$(basename "$sample" "_1.fq.gz")
    echo "Processing sample: ${base}"

    # Check if the paired reverse FASTQ file exists; skip if not
    if [ ! -f "${dirfastq}/${base}_2.fq.gz" ]; then
        echo "Error: Paired file not found: ${dirfastq}/${base}_2.fq.gz"
        continue
    fi

    # Initialize success flag for the sample processing
    success=true

    # Run the processing pipeline in a subshell and as a background job
    {
        # Alignment using Bowtie2
        echo "Running bowtie2 alignment..."
        bowtie2 -p "$BOWTIE_THREADS" -I 75 -X 1000 --no-discordant --no-mixed \
            -x "$ref_index" \
            -1 "${dirfastq}/${base}_1.fq.gz" \
            -2 "${dirfastq}/${base}_2.fq.gz" \
            -S "${dirsam}/${base}.sam" \
            2> "${dirsam}/${base}.bowtie2.v4.prefixchr.log"
        if [ $? -ne 0 ]; then
            echo "Bowtie2 alignment failed for sample: ${base}"
            success=false
        fi

        # Convert SAM to BAM (if previous step succeeded)
        if [ "$success" = true ]; then
            echo "Converting SAM to BAM..."
            samtools view -bS -h -@ "$SAMTOOLS_THREADS" "${dirsam}/${base}.sam" > "${dirsam}/${base}.bam"
            if [ $? -ne 0 ]; then
                echo "SAM to BAM conversion failed for sample: ${base}"
                success=false
            fi
        fi
      
        # Filter BAM by quality and flags (if previous step succeeded)
        if [ "$success" = true ]; then
            echo "Filtering BAM file..."
            samtools view -F "$FILTER_FLAG" -q "$QUALITY_THRESHOLD" -@ "$SAMTOOLS_THREADS" -b "${dirsam}/${base}.bam" > "${dirsam}/${base}.filter.bam"
            if [ $? -ne 0 ]; then
                echo "BAM filtering failed for sample: ${base}"
                success=false
            fi
        fi
      
        # Sort the filtered BAM (if previous step succeeded)
        if [ "$success" = true ]; then
            echo "Sorting BAM file..."
            samtools sort -@ "$SAMTOOLS_THREADS" "${dirsam}/${base}.filter.bam" -o "${dirsam}/${base}.filter.sorted.bam"
            if [ $? -ne 0 ]; then
                echo "BAM sorting failed for sample: ${base}"
                success=false
            fi
        fi
       
        # Mark and remove duplicates using Sambamba (if previous step succeeded)
        if [ "$success" = true ]; then
            echo "Marking duplicates..."
            sambamba markdup -r -t "$SAMBAMBA_THREADS" -p "${dirsam}/${base}.filter.sorted.bam" "${dirsam}/${base}.filter.sorted.dedup.bam"
            if [ $? -ne 0 ]; then
                echo "Duplicate marking failed for sample: ${base}"
                success=false
            fi
        fi

        # Generate flag statistics for the final BAM (if previous step succeeded)
        if [ "$success" = true ]; then
            echo "Generating BAM statistics..."
            samtools flagstat "${dirsam}/${base}.filter.sorted.dedup.bam" > "${dirsam}/${base}.bam-log"
            if [ $? -ne 0 ]; then
                echo "BAM statistics generation failed for sample: ${base}"
                success=false
            fi
        fi

        # If all steps succeeded, delete original FASTQ files (currently commented out)
        if [ "$success" = true ]; then
            echo "Deleting original FASTQ files..."
            #rm -f "${dirfastq}/${base}_1.fq.gz" "${dirfastq}/${base}_2.fq.gz"
        fi

        # Cleanup intermediate files regardless of success
        echo "Cleaning up intermediate files..."
        rm -f "${dirsam}/${base}.sam"
        rm -f "${dirsam}/${base}.bam"
        rm -f "${dirsam}/${base}.filter.bam"
        rm -f "${dirsam}/${base}.filter.sorted.bam"

        echo "Completed processing sample: ${base}"
    } &  # Run the subshell in the background

    # Increment the job counter
    job_count=$((job_count + 1))

    # If the maximum concurrent jobs is reached, wait for them to finish
    if [ "$job_count" -ge "$max_jobs" ]; then
        wait  # Wait for all background jobs to complete
        job_count=0  # Reset the job counter
    fi
done

# Wait for any remaining background jobs to complete
wait

# Final status message
echo "All samples processed successfully!"
