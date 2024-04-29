# MnM: A pan-cancer classifier for pediatric malignancies

## What is this?

`MnM` is an R-package that includes the code for running the newly developed machine-learning algorithm M&M, as short for Minority & Majority classifier. This classifier enables the classification of pediatric tumor (sub)types based on RNA-seq data. Besides the M&M classifier, `MnM` contains many functions to visualize M&M's performance on reference and test cohorts.

## Why and when to use this package?

With over 120 documented tumor types and 250 tumor subtypes, identifying the correct malignancy during the diagnostic procedure within a hospital remains a challenging but crucial process within pediatric oncology. M&M allows for inclusion of many rare pediatric malignancies with as few as three available samples, occurring with rates of less than once for every 500 children diagnosed with cancer. With the development of M&M, pathologists can be assisted during the diagnostic process to reduce inter-observer variability and help recognizing rare pediatric malignancies.

## Install

```{r}
library(remotes)
remotes::install_github("princessmaximacenter/MnM/", dependencies = T)

```

## Usage

Please see the supplied tutorial, vignettes, and documentation within R on how to properly use all functions. Pre-trained models can be obtained from (LOCATION).

## Contact

In case of questions, suggestions or additional comments, please reach out to us via (EMAIL?).
