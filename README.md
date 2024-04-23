# MnM: A pan-cancer classifier for pediatric malignancies

## What is this?

With over 120 documented tumor types and 250 tumor subtypes, identifying the correct malignancy during the diagnostic procedure within a hospital remains a challenging but crucial process within pediatric oncology.

`MnM` is an R-package that includes the code for running the newly developed machine-learning algorithm M&M, as short for Minority & Majority classifier. This classifier enables the classification of pediatric tumor (sub)types based on RNA-seq data. Within the setup, M&M allows for inclusion of many rare pediatric malignancies with as few as three available samples, to help with the tumor (sub)types diagnoses that are otherwise easily missed and prone to inter-observer variability. Besides the M&M classifier, `MnM` contains many functions to visualize M&M's performance.

## Install

```{r}
library(remotes)
remotes::install_github("princessmaximacenter/MnM/")

```

## Usage

Please see the supplied tutorial, vignettes, and documentation within R on how to properly use all functions. Pre-trained models can be obtained from (LOCATION).

## Contact

In case of questions, suggestions or additional comments, please reach out to us via (EMAIL?).
