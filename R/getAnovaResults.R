#' Calculate F-statistic ANOVA for feature selection
#'
#' This function calculates the differences between the RNA-transcripts expression of all tumor subtypes
#' from the training data using an analysis of variance (ANOVA) setup, assuming equal variance.
#' In this way, the RNA-transcripts that are most variable between
#' the tumor subtypes can be selected, based on their F-scores.
#'
#' @param allMaterial Dataframe containing both the cleaned RNA-seq count data and the metadata
#' containing the tumor subtype labels. This dataframe is a concatenation of the count-data and metadata files used before.
#' @param allGenes All RNA-transcripts that should be compared between tumor subtypes using the ANOVA-setup.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#'
#' @return Dataframe containing the genes with their associated F-scores and p-scores from ANOVA.
#'
getAnovaResults <- function(allMaterial,
                            allGenes,
                            classColumn) {
  anovaValues <- data.frame(allGenes = allGenes, p_val = NA, F_val = NA)
  for (i in seq(1:base::length(allGenes))) {
    ANOVARes <- stats::aov(allMaterial[,allGenes[i]] ~ allMaterial[, classColumn])
    anovaValues[i,"F_val"] <- summary(ANOVARes)[[1]][["F value"]][1]
    anovaValues[i,"p_val"] <- summary(ANOVARes)[[1]][["Pr(>F)"]][1]
  }
  return(anovaValues)
}
