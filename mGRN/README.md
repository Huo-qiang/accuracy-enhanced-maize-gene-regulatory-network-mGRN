# ampDAP-seq （DAP-seq） Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)  
*(简要中文描述：这是一个DAP-seq分析管道的Bash脚本集合，从reads比对到motif发现和比较。脚本已注解，适合SLURM集群运行。)*

## Overview

This repository contains a series of annotated Bash scripts that form a complete pipeline for analyzing DAP-seq (DNA Affinity Purification sequencing) data or similar high-throughput sequencing assays (e.g., ChIP-seq, ATAC-seq). The pipeline processes raw sequencing data through alignment, peak calling, quality control, sequence extraction, and motif analysis using tools like Bowtie, MACS3, Bedtools, phantompeakqualtools, MSPC, and the MEME suite.

The scripts are numbered sequentially (1 to 10) to indicate the recommended order of execution. Each script is self-contained, with detailed comments explaining its purpose, variables, and commands. They support concurrency where appropriate (e.g., for large datasets) and include error handling, output skipping (if results exist), and optional SLURM directives for cluster environments.

### Pipeline Flow
1. **Read Alignment**: Align reads to a reference genome.
2. **Peak Calling**: Identify peaks using MACS3.
3. **File Conversion**: Convert BAM to BigWig for visualization.
4. **Quality Metrics**: Calculate FRiP scores.
5. **QC Analysis**: Run phantompeakqualtools for BAM quality control.
6. **Multi-Sample Peak Calling**: Use MSPC for consensus peaks across replicates.
7. **Sequence Extraction**: Extract FASTA sequences from top peaks.
8. **Motif Discovery**: Run MEME-ChIP to discover de novo motifs.
9. **Motif Enrichment**: Perform enrichment analysis with AME.
10. **Motif Comparison**: Compare motifs to known databases using TomTom.

## Requirements

- **Operating System**: Linux (tested on Ubuntu or similar; SLURM support for clusters).
- **Tools** (install via conda/mamba, apt, or manually):
  - Bowtie2 (for alignment).
  - MACS3 (for peak calling).
  - Samtools and Bedtools (for BAM/BED manipulation and FASTA extraction).
  - phantompeakqualtools (for QC; requires Rscript).
  - MSPC (Multi-Sample Peak Caller).
  - MEME suite (version 5.5.1 or later; for MEME-ChIP, AME, TomTom).
  - Other: awk, samtools, bedtools, R (for phantompeakqualtools).
- **Reference Files**: Genome FASTA (e.g., Zea_mays.B73_RefGen_v4.49.dna.chr.fa), motif databases (e.g., JASPAR2024, ArabidopsisDAPv1).
- **Environment**: SLURM cluster recommended for large datasets; adjust paths and variables in scripts as needed.
- **Installation Tip**: Use conda/mamba to create an environment: `mamba create -n dap_pipeline bowtie2 macs3 bedtools samtools r-base meme`.

### Scripts Overview

| Script                  | Description |
|-------------------------|-------------|
| **1-bowtie.sh**        | Aligns FASTQ reads to a reference genome using Bowtie2, with deduplication and sorting. Outputs BAM files. |
| **2-macs3.sh**         | Performs peak calling on BAM files using MACS3, handling replicates and controls. Outputs BED files. |
| **3-bam2bw.sh**        | Converts BAM files to BigWig format for genome browser visualization using deepTools or similar. |
| **4-FRiP.sh**          | Calculates Fraction of Reads in Peaks (FRiP) scores using Samtools and Bedtools. |
| **5-Phantompeakqualtools_QC.sh** | Runs phantompeakqualtools for BAM quality control, generating PDF plots and text reports. Supports concurrency. |
| **6-MSPC.sh**          | Runs MSPC for multi-sample consensus peak calling on BED file pairs (e.g., replicates). Outputs consensus BED. |
| **7-bedtools_fasta.sh**| Extracts FASTA sequences from top peaks (e.g., top 600, centered and extended) using Bedtools. Supports high concurrency. |
| **8-memechip.sh**      | Discovers de novo motifs in FASTA files using MEME-ChIP, with parameters for width and skipping analyses. |
| **9-AME.sh**           | Performs motif enrichment analysis on FASTA files using AME against two databases. Supports concurrency. |
| **10-TOMTOM.sh**       | Compares discovered motifs to known databases using TomTom, identifying similarities. |

### Example Workflow
- Input: Raw FASTQ files.
- Output: Aligned BAMs → Peaks (BED) → QC metrics → Extracted FASTA → Motifs (MEME) → Enrichment/Comparison results.
- Total runtime depends on dataset size; use SLURM for parallel processing.

## Troubleshooting
- **Errors**: Check paths, tool installations, and permissions. Scripts include error handling (e.g., skipping missing files).
- **Customization**: Adjust parameters like `max_jobs` for concurrency or thresholds in peak calling/motif tools.
- **Logs**: Each script echoes progress and errors to stdout.

## Contributing
Pull requests are welcome! For major changes, open an issue first to discuss.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. (Create a LICENSE file if needed.)

## Acknowledgments
- Built with tools from the MEME suite, MACS3, Bedtools, and others.
- Inspired by standard DAP-seq workflows for plant genomics (e.g., maize/Arabidopsis).
