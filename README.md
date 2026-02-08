## Overview

Cancer neuroscience highlights the critical role of neural signaling in tumor biology, yet a comprehensive pan-cancer understanding of neuroregulatory dysregulation remains lacking. Here, we systematically characterized 130 neurotransmitter receptor (NTR) genes across 9,125 tumors from 33 cancer types, and performed validation using multiple independent cohorts. This repository aims to share the analytical workflows, processed data, and supporting resources to facilitate deeper understanding of this study, enable reproducibility, and support reuse of results of interest by readers and reviewers. Detailed package version specifications, complete parameter settings for the major analytical procedures, and guidance for the core analysis steps are provided in the README.md file. Additional scripts implementing specific functions can be accessed by referencing the corresponding script files.

## Repo Contents

* [code](https://github.com/BioinfoDriver/CancerNeuro/tree/main/code): tidy R functions and R script.
* [data](https://github.com/BioinfoDriver/CancerNeuro/tree/main/data): original and preprocessed data used for analysis and share.
* [result](https://github.com/BioinfoDriver/CancerNeuro/tree/main/result): important middle results and final results of manuscript.

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
