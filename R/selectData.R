#' Title
#'
#' @param refCohort
#' @param select_vector
#' @param depth
#' @param classColumn
#' @param n
#'
#' @return
#' @export
#'
#' @examples
selectData <- function(refCohort, select_vector = NA, depth = NA, classColumn, n = 0) {
  if (is.na(select_vector) & is.na(depth)) {
    metaData <- refCohort$metaData
    colnumDepth <- which(colnames(refCohort$metaData) == classColumn)
    abundantTumorTypes <- names(table(metaData[,colnumDepth]))[table(metaData[,colnumDepth]) > n]
    metaData <- metaData[metaData[,colnumDepth] %in% abundantTumorTypes,]
  } else {
    metaData <- refCohort$metaData[refCohort$metaData[, depth] %in% select_vector,]
    colnumDepth <- which(colnames(refCohort$metaData) == depth)

    abundantTumorTypes <- names(table(metaData[,colnumDepth+1]))[table(metaData[,colnumDepth+1]) > n]
    metaData <- metaData[metaData[,colnumDepth+1] %in% abundantTumorTypes,]
  }
  all_RNA_seq_biomaterial_IDs <- refCohort$counts %>%
    colnames()

  selected_RNA_seq_biomaterial_IDs <- rownames(metaData)

  countData <- refCohort$counts[, all_RNA_seq_biomaterial_IDs %in% selected_RNA_seq_biomaterial_IDs]

  refCohort$counts <- countData
  refCohort$metaData <- metaData

  return(refCohort)
}
