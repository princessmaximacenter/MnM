#' Clean reference cohort
#'
#' Function to remove duplicate samples from the same patient from the reference cohort,
#' and tumor (sub)types for which there are not enough samples available (< _n_).
#'
#' @param refCohort R-object containing both the RNA-seq count data matrix ($counts) and the associated metadata ($metaData).
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param n What is the minimum number of samples that is required per tumor (sub)type?
#' @param inputDir Directory where the CSR-data is stored that is needed to check whether
#' there is only one data entry per patient per tumor type.
#' @export
#' @return R-object containing the cleaned RNA-seq count data matrix ($counts) and the associated metadata ($metaData).
#'
cleanData <- function(refCohort, classColumn = "Disease_sub_specification1", n = 3,
                      inputDir = "~/surfdrive/Shared/Kemmeren group/Research_Projects/RNA_classification_FW/data/input/CSR/2023_03/") {
  newRefCohort <- removeDuplicates(refCohort = refCohort,
                                   inputDir = inputDir)
  newRefCohort <- selectData(newRefCohort,
                             classColumn = classColumn,
                             n = n)
  return(newRefCohort)
}
