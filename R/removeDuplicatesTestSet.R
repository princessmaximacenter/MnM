#' Title
#'
#' @param newRefCohort
#' @param testSet
#'
#' @return
#' @export
#'
#' @examples
removeDuplicatesTestSet <- function(newRefCohort,
                                    testSet
) {
  patientsRefCohort <- unique(newRefCohort$metaData$PMCID)
  testSetMetaData <- testSet$metaData %>% filter(PMCID %notin% patientsRefCohort,
                                                 rownames(.) %notin% rownames(newRefCohort$metaData))

  all_RNA_seq_biomaterial_IDs <- testSet$counts %>%
    colnames()
  selected_RNA_seq_biomaterial_IDs <- rownames(testSetMetaData)
  testCountData <- testSet$counts[, all_RNA_seq_biomaterial_IDs %in% selected_RNA_seq_biomaterial_IDs]

  testSet$counts <- testCountData
  testSet$metaData <- testSetMetaData
  return(testSet)
}
