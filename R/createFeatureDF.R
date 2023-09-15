#' Extract importance scores from RF-models
#'
#' This function is designed to extract all the mean decreases in accuracy
#' from a performed weighted Random Forest (RF). The mean decreases in accuracy / Gini score
#' from all different #' genes selected in the ANOVA-procedure is calculated,
#' for all models within the modelList.
#'
#' @param modelList List containing different weighted RF models with importance scores.
#' @param whichAccuracyMeasure Which statistic calculated in the importance score of the RF
#' models do you want to use, the mean decrease in accuracy (MeanDecreaseAccuracy)
#' or the mean decrease in the Gini score (MeanDecreaseGini)?
#' @param nANOVAgenes How many genes did we select during the ANOVA procedure?
#'
#' @return Dataframe containing the mean decreases in accuracy / Gini score from removing
#' different genes during a weighted RF procedure for different models.
#' The genes are in the rows, the different models in the columns.
#' @export
#'
createFeatureDF <- function(modelList, whichAccuracyMeasure, nANOVAgenes) {
  ##Create an empty dataframe
  empty_df <- data.frame()

  #Recursively extract column values from ModelList
  for (i in seq_along(modelList)) {

    columnValues <- as.data.frame(modelList[[i]][["importance"]]) %>%
      select(all_of(whichAccuracyMeasure)) %>%
      arrange(desc(all_of(whichAccuracyMeasure))) %>%
      slice(1:nANOVAgenes) %>%
      rownames_to_column(var = "Gene") %>%
      mutate(model = i) %>%
      rename(value = whichAccuracyMeasure)

    #Append to Dataframe
    empty_df <- rbind(empty_df, columnValues)
  }

  #Create create a wide DF for comparing features
  accuracyValuesDF <- empty_df %>% pivot_wider(names_from = "model")

  #Store columnName as rownames
  accuracyValuesDF <- accuracyValuesDF %>% remove_rownames %>% column_to_rownames(var="Gene")

  return(accuracyValuesDF)
}
