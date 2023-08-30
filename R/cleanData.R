#' Title
#'
#' @param refCohort
#' @param classColumn
#' @param n
#' @param inputDir
#'
#' @return
#' @export
#'
#' @examples
cleanData <- function(refCohort, classColumn = "Disease_sub_specification1", n = 2,
                      inputDir = "~/surfdrive/Shared/Kemmeren group/Research_Projects/RNA_classification_FW/data/input/CSR/2023_03/") {
  newRefCohort <- removeDuplicates(refCohort = refCohort,
                                   inputDir = inputDir)
  newRefCohort <- selectData(newRefCohort,
                             classColumn = classColumn,
                             n = n)
  newRefCohort <- removeNotMalignantNeuroAndInflammation(newRefCohort = newRefCohort)
  return(newRefCohort)
}
