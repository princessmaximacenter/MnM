#' Remove bias from test set
#'
#' Remove samples belonging to patients that are present in the reference cohort from the test set,
#' and are thereby violating the assumption of independence between the reference cohort and test set.
#'
#' @param newRefCohort Reference cohort R-object containing both the RNA-seq count data (refCohort$counts) and
#' the metadata associated to this count data (refCohort$metaData).
#' @param testSet Testset R-object containing both the RNA-seq count data (testSet$counts) and
#' the metadata associated to this count data (testSet$metaData).
#'
#' @return The test set without samples belonging to patients from the reference cohort.

#' @export
#'
removeBiasTestSet <- function(newRefCohort,
                                    testSet
) {
  patientsRefCohort <- unique(newRefCohort$metaData$PMCID)
  testSetMetaData <- testSet$metaData %>% filter(PMCID %notin% patientsRefCohort,
                                                 rownames(.) %notin% rownames(newRefCohort$metaData))

  allRNAseqBiomaterialIDs <- testSet$counts %>%
    colnames()
  selectedRNAseqBiomaterialIDs <- rownames(testSetMetaData)
  testCountData <- testSet$counts[, allRNAseqBiomaterialIDs %in% selectedRNAseqBiomaterialIDs]

  testSet$counts <- testCountData
  testSet$metaData <- testSetMetaData
  return(testSet)
}
