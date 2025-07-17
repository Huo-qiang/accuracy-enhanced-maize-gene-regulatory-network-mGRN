# mGRN+ - Enhanced Maize Gene Regulatory Network

## Overview

This repository contains data, code, and analyses associated with our large-scale profiling of maize transcription factor (TF) binding sites and the construction of an accuracy-enhanced maize gene regulatory network (mGRN+). The network integrates chromatin accessibility and gene expression data to improve predictions of gene functions and key regulators in maize.

## Research Summary

Understanding gene regulatory networks (GRNs) is essential for improving maize yield and quality through molecular breeding approaches. The lack of comprehensive transcription factor (TF)-DNA interaction data has hindered accurate GRN predictions, limiting our insight into the regulatory mechanisms. 

In this work:

- We performed large-scale profiling of maize TF binding sites
- We obtained and collected reliable binding profiles for 513 TFs
- We identified 394,136 binding sites
- We constructed an accuracy-enhanced maize GRN (mGRN+) by integrating chromatin accessibility and gene expression data
- The mGRN+ comprises 397,699 regulatory relationships
- We divided the mGRN+ into multiple modules across six major tissues
- Using machine learning algorithms, we optimized the mGRN+ network to improve prediction accuracy of gene functions and key regulators
- Through independent genetic validation experiments, we confirmed the reliability of these predictions

## Repository Structure

- **ampDAP-seq/**: Contains scripts and data related to DNA affinity purification sequencing analyses
- **Motif/**: Contains MEME files for each TF's binding motif
- **TF_peaks_files/**: Contains peak files for each transcription factor
- **mGRN/**: Contains specific regulatory relationships and tissue-specific mGRN regulatory relationships across 6 tissues
- **mGRN construct/**: Scripts and methods used to construct the network
- **machine learning/**: Implementation of machine learning algorithms used to predict key TFs
- **ck_puc57_greylist_merge_sorted_final.bed**: DAP-seq greylist file for filtering false positive peaks

## Data Availability

This repository provides:

1. The largest collection of experimental TF binding sites in maize
2. A highly optimized regulatory network (mGRN+)
3. Tissue-specific network modules across 6 major tissues
4. Computational methods for network construction and optimization
5. TF binding motifs in MEME format

## Usage

### Prerequisites

- R (version ≥ 4.0)
- Python (version ≥ 3.6)
- Bioinformatics tools: Bowtie2, MACS3, MEME Suite, etc.

### Analysis Pipeline

The complete analysis pipeline is available in the repository, organized in sequential Bash scripts:

1. Read alignment with Bowtie2
2. Peak calling with MACS3
3. BAM to BigWig conversion for visualization
4. Quality metrics calculation (FRiP scores)
5. Quality control with phantompeakqualtools
6. Multi-sample peak calling with MSPC
7. Sequence extraction from top peaks
8. Motif discovery with MEME-ChIP
9. Motif enrichment analysis with AME
10. Motif comparison with TomTom

### Network Construction

To construct the mGRN+ network:

1. Follow the scripts in `mGRN construct/`
2. Integrate chromatin accessibility data
3. Combine with gene expression data
4. Apply the machine learning models in `machine learning/` for key TF prediction

## Applications

This work provides valuable resources for:
- Study of maize gene function
- Identification of key regulatory genes
- Crop improvement applications
- Molecular breeding approaches
- Investigation of tissue-specific regulatory mechanisms

## Citation

If you use the data or methods from this repository, please cite our paper:
Exploring maize transcriptional regulatory landscape through large-scale profiling of transcription factor binding sites.

## Contact

For questions or collaborations, please contact Qiang Huo [huoqiang@cau.edu.cn].

*This repository contains the data and methods supporting our research on maize gene regulatory networks, providing a comprehensive resource for the scientific community working on crop improvement and plant regulatory genomics.*
