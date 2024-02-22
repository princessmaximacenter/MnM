#' Obtain empty tiles for confusion matrix
#'
#' Get tiles from the confusion matrix that are never used for a prediction.
#' These tiles are important for the generation of the eventual confusion matrix.
#'
#' @param predictionsMM Predictions for the tumor (sub)types with their associated probability scores.
#' @param probabilityScore Probability score that you would like to use as a threshold for calling samples
#' 'classifiable' or 'non-classifiable'.
#' @param abbreviations Dataframe containing the abbreviations and full names of the tumor (sub)types, for the eventual confusion matrix plot labels.
#' The abbreviation should be stored in a column with header 'abbreviation'.
#' @param classColumn Specifies the name of the column of the abbreviations dataframe that contains the full tumor type label.
#'
#' @return A dataframe containing the empty tiles for a confusion matrix plot.
#' @export
#' @import caret
getNonAvailableTiles <- function(predictionsMM,
         probabilityScore,
         abbreviations,
         classColumn) {


  predictionsMMFiltered <- predictionsMM %>% filter(probability1 > probabilityScore)

  tumorConfusionMatrix <- confusionMatrix(factor(predictionsMMFiltered$predict,
                                                 levels = unique(c(predictionsMM$originalCall, predictionsMM$predict))),
                                          factor(predictionsMMFiltered$originalCall, levels = unique(c(predictionsMM$originalCall, predictionsMM$predict))),
                                          dnn = c("Prediction", "Reference"))
  predictionFrequencies <- tumorConfusionMatrix$table %>% as.data.frame()
  nonAvailableTiles <- predictionFrequencies %>% filter(Freq == 0, Prediction != Reference)

  nonAvailableTiles$Domain <- NA
  nonAvailableTiles$Prediction <- as.character(nonAvailableTiles$Prediction)
  nonAvailableTiles$Reference <- as.character(nonAvailableTiles$Reference)

  for (i in seq(1:nrow(nonAvailableTiles))) {

    if (nonAvailableTiles$Prediction[i] %in% abbreviations[, classColumn]) {
      nonAvailableTiles$Prediction[i] <- abbreviations[abbreviations[,classColumn] == nonAvailableTiles$Prediction[i], "abbreviation"]
    }
    if (nonAvailableTiles$Reference[i] %in% abbreviations[, classColumn]) {
      nonAvailableTiles$Reference[i] <- abbreviations[abbreviations[, classColumn] == nonAvailableTiles$Reference[i], "abbreviation"]
    }
  }
  return(nonAvailableTiles)

}
