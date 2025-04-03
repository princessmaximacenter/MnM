# M&M: A pan-cancer classifier for paediatric tumours

## What is MnM?

`MnM` is an R-package that includes the code for running the newly developed machine-learning algorithm M&M, as short for Minority & Majority classifier. This classifier enables the classification of paediatric tumour (sub)types based on RNA-seq data. Besides the creation of the classifier models, `MnM` contains many functions to visualize M&M's performance on reference and test cohorts.

## Why and when to use this package?

With over 120 documented tumour types and 250 tumour subtypes, identifying the correct tumour entity during the diagnostic procedure remains a challenging but crucial process within paediatric oncology. M&M allows for inclusion of many rare paediatric tumours with as few as three available samples, occurring with rates of less than once for every 500 children diagnosed with cancer. M&M is designed to assist pathologists during the diagnostic process to reduce inter-observer variability and help recognising rare paediatric tumours.

## Install

```{r}
library(remotes)
remotes::install_github("princessmaximacenter/MnM", dependencies = T)

```

## Usage

Please see the supplied tutorial, and documentation within R on how to properly use all functions. If you wish to familiarize yourself with `MnM`, the tutorial is the best starting point (DataTutorial). For this purpose, you can download the ZIP of this github page (click <> Code, Download ZIP). 

RNA TPM-normalized count data and metadata of the reference cohort and test cohort can be obtained from ArrayExpress (accession E-MTAB-14038). The final resulting pre-trained models, which can be used to classify new incoming RNA-samples, can be obtained from Zenodo (https://zenodo.org/records/14167359). These models now are also capable of performing missing gene imputation, making them more user-friendly. Please note that RNA-transcript count rownames are required to be HGNC-symbols, ENSEMBL-IDs currently cannot be used. 

## Contact

In case of questions, suggestions or additional comments, please reach out to us via **p.kemmeren[AT]prinsesmaximacentrum.nl** (orcid ID: 0000-0003-2237-7354).
