#' Calculate F-statistic ANOVA for gene selection
#'
#' This function calculates the differences between the gene expression of all tumor (sub)types
#' from the reference cohort using an ANOVA setup. In this way, the genes that are most variable between
#' the tumor (sub)types can be selected, based on their F-scores.
#'
#' @param allMaterial Dataframe containing both the cleaned RNA-seq count data and the metadata
#' containing the tumor (sub)type labels. This dataframe is a concatenation of the count-data and metadata files used before.
#' @param allGenes All genes that should be compared between tumor (sub)types using the ANOVA-setup.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#'
#' @return Dataframe containing the genes with their associated F-scores and p-scores from ANOVA.
#' @export
#'
#' @examples
getAnovaResults <- function(allMaterial, allGenes, classColumn) {
  anovaValues <- data.frame(allGenes = allGenes, p_val = NA, F_val = NA)
  for (i in seq(1:length(allGenes))) {
    ANOVARes <- aov(allMaterial[,allGenes[i]] ~ allMaterial[, classColumn])
    anovaValues[i,"F_val"] <- summary(ANOVARes)[[1]][["F value"]][1]
    anovaValues[i,"p_val"] <- summary(ANOVARes)[[1]][["Pr(>F)"]][1]
  }
  return(anovaValues)
}
