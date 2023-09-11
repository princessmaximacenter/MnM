#' Selecting tumor (sub)types
#'
#'This function removes all tumor (sub)types from the dataset which do not have _n_ or more samples. The M&M setup
#'requires an absolute minimum of 3 samples per tumor type for the cross-validation setup,
#' and 2 for the setup for new samples due to the use of the F-statistic for gene selection.
#' However, larger minimum numbers can be specified within this function as well.
#'
#' @param refCohort Reference cohort R-object containing both the RNA-seq count data (refCohort$counts) and
#' the metadata associated to this count data (refCohort$metaData).
#' @param classColumn Column in the metadata where the prediction labels for the classification process are specified.
#' @param n Specification of the minimum number of samples that should be present for each tumor (sub)type.
#'
#' @return The new reference cohort containing only the tumor (sub)types with more than _n_ entries within the dataset.
#' @export
#'
#'@import tidyverse dplyr magrittr
#'
selectData <- function(refCohort,
                       classColumn,
                       n = 3) {

  metaData <- refCohort$metaData
  colnumDepth <- which(colnames(refCohort$metaData) == classColumn)
  abundantTumorTypes <- names(table(metaData[,colnumDepth]))[table(metaData[,colnumDepth]) >= n]
  metaData <- metaData[metaData[,colnumDepth] %in% abundantTumorTypes,]

  allRNAseqBiomaterialIDs <- refCohort$counts %>%
    colnames()

  selectedRNAseqBiomaterialIDs <- rownames(metaData)

  countData <- refCohort$counts[, allRNAseqBiomaterialIDs %in% selectedRNAseqBiomaterialIDs]

  refCohort$counts <- countData
  refCohort$metaData <- metaData

  return(refCohort)
}
