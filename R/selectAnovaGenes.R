#' Select most informative genes Minority Classifier
#'
#' This function allows for the extraction of the genes that are most differently expressed between the different tumor (sub)types.
#' The gene selection is based on an analysis of variance (ANOVA) with equal variances, where the genes with the highest F-scores are selected.
#'
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis for the reference cohort.
#' @param countDataRef Matrix containing the RNA-transcript per million data. Patients are in the columns,
#' different genes in the rows.
#' @param nANOVAgenes How many genes should we select using the F-statistic of ANOVA?
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#'
#' @return Vector containing the names of the most interesting genes according to the F-statistic of ANOVA.
#'
selectAnovaGenes <- function(metaDataRef,
                             countDataRef,
                             nANOVAgenes,
                          classColumn
) {

  countDataRef <- t(countDataRef) %>% as.data.frame
  allGenes <- colnames(countDataRef)

  classesWith2 <- table(metaDataRef[,classColumn])[table(metaDataRef[,classColumn]) == 2] %>%
    names(.)

  if (length(classesWith2) > 0 ) {
    countDataRef$class <- as.factor(metaDataRef[rownames(countDataRef),classColumn])
    countDataRef <- createExtraData(countDataRef, classesWith2)
    metaDataRef$class <- as.character(metaDataRef[,classColumn])
    metaDataRef <- createExtraMetaData(metaDataRef = metaDataRef, classesWith2 = classesWith2)
  }

  countDataRef$Patient <- rownames(countDataRef)

  metaDataRef$Patient <- rownames(metaDataRef)
  allMaterial <- left_join(metaDataRef, countDataRef, by = "Patient")
 # allGenes <- allMaterial %>%
  #  select(where(is.numeric)) %>% colnames

  # Perform an ANOVA test on all genes
  results <- getAnovaResults(allMaterial, allGenes, classColumn)

  # Arrange the genes so that the genes with the highest F-scores are at the top
  filterResults <- results %>% arrange(desc(F_val))

  # Select the top n F-score genes
  interestingAnovaGenes <- head(filterResults$allGenes, n = nANOVAgenes)
  return(interestingAnovaGenes)
}
