## Overview

Cancer neuroscience highlights the critical role of neural signaling in tumor biology, yet a comprehensive pan-cancer understanding of neuroregulatory dysregulation remains lacking. Here, we systematically characterized 130 neurotransmitter receptor (NTR) genes across 9,125 tumors from 33 cancer types, and performed validation using multiple independent cohorts. This repository aims to share the analytical workflows, processed data, and supporting resources to facilitate deeper understanding of this study, enable reproducibility, and support reuse of results of interest by readers and reviewers. Detailed package version specifications, complete parameter settings for the major analytical procedures, and guidance for the core analysis steps are provided in the README.md file. Additional scripts implementing specific functions can be accessed by referencing the corresponding script files.

## Repo Contents

* [code](https://github.com/BioinfoDriver/CancerNeuro/tree/main/code): This directory contains the organized R functions and scripts used throughout the study, structured into subdirectories corresponding to each main section of the manuscript:
   - SomaticMutationLandscape — corresponds to Section 1
   - SomaticCopyNumberAlterations & ConditionalSelectionDependency — correspond to Section 2
   - TranscriptomicDysregulation & PrognosticSignificance — correspond to Section 3
   - EpigeneticRegulation — corresponds to Section 4
   - NeuroregulatorySubtypes — corresponds to Section 5
   - Subtype4LggSubtype & Subtype1LihcSubtype — correspond to Section 6
  - Each subdirectory includes scripts to reproduce the analyses and figures presented in the respective section.
* [data](https://github.com/BioinfoDriver/CancerNeuro/tree/main/data): This directory stores the original and preprocessed data used in the analyses. Due to file size limitations, some datasets are provided via external download links. If you require access to data not directly available here, please contact the lead author Yajing Zhang (zhangyajing@hrbmu.edu.cn) with a reasonable request.
* [result](https://github.com/BioinfoDriver/CancerNeuro/tree/main/result): This directory contains key intermediate and final outputs from the analyses, including processed tables, statistical summaries, and visualization-ready files that support the findings presented in the manuscript.

## Instructions for Core Analytical Procedures
### Annotating Genomic Point Mutations and Short Indels
* Requirements
   - Java: Version 1.8 
   - GATK4: 4.6.2.0

* Execution Command
   - ./gatk Funcotator \
   - --variant variants.vcf \
   - --reference Homo_sapiens_assembly19.fasta \
   - --ref-version hg19 \
   - --data-sources-path funcotator_dataSources.v1.8.hg19.20230908s\
   - --output variants.funcotated.maf \
   - --output-file-format MAF
   - --ignore-filtered-variants false \
   - --transcript-selection-mode CANONICAL

* Parameter Description
   - --variant: Path to the input VCF file containing variants to be annotated.
   - --reference: Path to the reference sequence (FASTA format) used for alignment.
   - --ref-version: Human reference genome version (e.g. hg19, hg38, etc.).
   - --data-sources-path: Path to the Funcotator data sources directory. Data can be downloaded from [Broad Institute Google Cloud] ([https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator).
   - --output: Output file path for annotated variants.
   - --output-file-format: Output format: VCF, MAF, or SEG.
   - --ignore-filtered-variants: Ignore/drop variants that have been filtered in the input. These variants will not appear in the output file.
   - --transcript-selection-mode: Method of detailed transcript selection. This will select the transcript for detailed annotation (CANONICAL, ALL, or BEST_EFFECT).
   - For detailed usage and tutorials, please refer to the official GATK documentation:
[Funcotator Information and Tutorial](
https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial).

## Instructions for Use
For readers who want to obtain raw/result data, locate data file, then download it with one of following ways:

* In Github, download file by clicking either `Download` button or `Raw` button at corresponding data page

* Use linux command `wget` or `curl`, fo example, you can download neurotransmitter receptors by

  wget `https://github.com/BioinfoDriver/CancerNeuro/tree/main/data/neurotransmitterReceptors.rds`

Or you can download whole respository with one of following ways:

* Clone this repository with `git clone https://github.com/BioinfoDriver/CancerNeuro.git`

* Download whole respository by clicking `Download` button at top right of url page https://github.com/BioinfoDriver/CancerNeuro

## Reproduce analysis

For readers who want to reproduce analysis shown in manuscript, please [install R](https://cran.r-project.org/) in your computer.

## Test Environment
* System: **Linux**
* Software: **R v4.0.2**
