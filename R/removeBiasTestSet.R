#' Remove bias from test set
#'
#' Remove samples belonging to patients that are present in the reference cohort from the test set,
#' and are thereby violating the assumption of independence between the reference cohort and test set.
#' Also, the samples that can never be predicted correctly (labels not within reference cohort) are removed.
#'
#' @param newRefCohort Reference cohort R-object containing both the RNA-seq count data (refCohort$counts) and
#' the metadata associated to this count data (refCohort$metaData).
#' @param testSet Testset R-object containing both the RNA-seq count data (testSet$counts) and
#' the metadata associated to this count data (testSet$metaData).
#' @param classColumn Which column in the metadata file should be used to filter out data points
#' from the test set that do not have a corresponding label?
#' @return The test set without samples belonging to patients from the reference cohort.

#'
removeBiasTestSet <- function(metaDataRef,
                                    metaDataTest,
                              classColumn,
                              patientIDColumn
) {
  `%notin%` <- Negate(`%in%`)
  patientsRefCohort <- unique(metaDataRef[, patientIDColumn])
  metaDataTestNew <- metaDataTest%>% filter(!!sym(patientIDColumn) %notin% patientsRefCohort,
                                                 rownames(.) %notin% rownames(metaDataRef)
                                                 #!!sym(classColumn) %in% unique(newRefCohort$metaData[,classColumn])
                                                 )

  # allRNAseqBiomaterialIDs <- testSet$counts %>%
  #   colnames()
  # selectedRNAseqBiomaterialIDs <- rownames(testSetMetaData)
  # testCountData <- testSet$counts[, allRNAseqBiomaterialIDs %in% selectedRNAseqBiomaterialIDs]
  #
  # testSet$counts <- testCountData
  # testSet$metaData <- testSetMetaData


  return(metaDataTestNew)
}
