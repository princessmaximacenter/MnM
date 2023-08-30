#' Title
#'
#' @param newRefCohort reference cohort, cleaned for everything except the not malignant neuro and inflammation samples
#'
#' @return
#' @export newRefCohort reference cohort without the not malignant neuro and inflammation samples
#'
#' @examples
removeNotMalignantNeuroAndInflammation <- function(newRefCohort) {

  newRefCohort$metaData <- newRefCohort$metaData %>% filter(Disease_sub_class %notin% c("Not malignant neuro", "Not malignant inflammation"))
  all_RNA_seq_biomaterial_IDs <- newRefCohort$counts %>%
    colnames()

  selected_RNA_seq_biomaterial_IDs <- rownames(newRefCohort$metaData)

  newRefCohort$counts <- newRefCohort$counts[, all_RNA_seq_biomaterial_IDs %in% selected_RNA_seq_biomaterial_IDs]

  return(newRefCohort)
}
