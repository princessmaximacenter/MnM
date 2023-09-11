#' Title
#'
#' @param metaData
#' @param countData
#' @param nANOVAgenes
#' @param classColumn
#'
#' @return
#' @export
#'
getAnovaGenes <- function(metaData, countData, nANOVAgenes = 1000,
                          classColumn
) {

  countData <- t(countData) %>% as.data.frame
  classesWith2 <- table(metaData[,classColumn])[table(metaData[,classColumn]) == 2] %>%
    names(.)

  if (length(classesWith2) > 0 ) {
    countData$class <- as.factor(metaData[,classColumn])
    #countData$class <- as.factor(countData$class)
    countData <- createExtraData(countData, classesWith2)
    metaData$class <- as.character(metaData[,classColumn])
    metaData <- createExtraMetaData(metaData = metaData, classesWith2 = classesWith2)
  }

  countData$Patient <- rownames(countData)

  metaData$Patient <- rownames(metaData)
  allMaterial <- left_join(metaData, countData, by = "Patient")
  allGenes <- allMaterial %>%
    select(where(is.numeric)) %>% colnames

  # Perform an ANOVA test on all genes
  results <- getAnovaResults(allMaterial, allGenes, classColumn)

  # Arrange the genes so that the genes with the highest F-scores are at the top
  filterResults <- results %>% arrange(desc(F_val))

  # Select the top n F-score genes
  interestingAnovaGenes <- head(filterResults$allGenes, n = nANOVAgenes)
  return(interestingAnovaGenes)
}
