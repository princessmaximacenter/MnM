#' Select most informative RNA-transcripts as features for Minority classifier
#'
#' This function allows for the extraction of the RNA-transcripts that are most differently expressed between the different tumor subtypes.
#' The RNA-transcripts selection is based on an analysis of variance (ANOVA) with equal variances,
#' where the RNA-transcripts with the highest F-scores are selected.
#'
#' @param metaDataRef Metadata file containing the links between the samples and the tumor subtype diagnosis for the reference cohort.
#' @param countDataRef Matrix containing the RNA-transcript per million data. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param nANOVAgenes How many RNA-transcripts should we select using the F-statistic of ANOVA?
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#'
#' @return Vector containing the names of the most interesting RNA-transcripts according to the F-statistic of ANOVA.
#' @import utils
#'
selectAnovaGenes <- function(metaDataRef,
                             countDataRef,
                             nANOVAgenes,
                          classColumn
) {

  countDataRef <- base::t(countDataRef) %>% base::as.data.frame()
  allGenes <- base::colnames(countDataRef)

  classesWith2 <- base::table(metaDataRef[,classColumn])[base::table(metaDataRef[,classColumn]) == 2] %>%
    base::names(.)

  if (base::length(classesWith2) > 0 ) {
    countDataRef$class <- base::as.factor(metaDataRef[rownames(countDataRef),classColumn])
    countDataRef <- createExtraData(countDataRef, classesWith2)
    metaDataRef$class <- base::as.character(metaDataRef[,classColumn])
    metaDataRef <- createExtraMetaData(metaDataRef = metaDataRef,
                                       classesWith2 = classesWith2)
  }

  countDataRef$Sample <- base::rownames(countDataRef)

  metaDataRef$Sample <- base::rownames(metaDataRef)
  allMaterial <- dplyr::left_join(metaDataRef, countDataRef, by = "Sample")

  # Perform an ANOVA test on all RNA-transcripts
  results <- getAnovaResults(allMaterial = allMaterial,
                             allGenes = allGenes,
                             classColumn = classColumn)

  # Arrange the RNA-transcripts so that the RNA-transcripts with the highest F-scores are at the top
  filterResults <- results %>% dplyr::arrange(dplyr::desc(F_val))

  # Select the top n F-score RNA-transcripts
  interestingAnovaGenes <- utils::head(filterResults$allGenes,
                                n = nANOVAgenes)
  return(interestingAnovaGenes)
}
