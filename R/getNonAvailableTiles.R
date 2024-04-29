#' Obtain empty tiles for confusion matrix
#'
#' Get tiles from the confusion matrix that are never used for a classification.
#' These tiles are important for the generation of the eventual confusion matrix.
#'
#' @param predictionsMM Predictions for the tumor (sub)types with their associated probability scores.
#' @param probabilityThreshold Probability score that you would like to use as a threshold for calling samples
#' 'confident' or not.
#' @param abbreviations  Dataframe containing the links between the tumor subtype and their abbreviations ($abbreviationSubtype),
#' the tumor types and their abbreviations ($abbreviationTumorType), and the domain.
#' @param classColumn Specifies the name of the column of the abbreviations dataframe that contains the full tumor type label.
#'
#' @return A dataframe containing the empty tiles for a confusion matrix plot.
#' @import caret
getNonAvailableTiles <- function(predictionsMM,
         probabilityThreshold,
         abbreviations,
         classColumn) {


  predictionsMMFiltered <- predictionsMM %>% dplyr::filter(probability1 > probabilityThreshold)

  tumorConfusionMatrix <- caret::confusionMatrix(factor(predictionsMMFiltered$predict,
                                                 levels = unique(c(predictionsMM$originalCall, predictionsMM$predict))),
                                          factor(predictionsMMFiltered$originalCall, levels = unique(c(predictionsMM$originalCall, predictionsMM$predict))),
                                          dnn = c("Prediction", "Reference"))
  predictionFrequencies <- tumorConfusionMatrix$table %>% as.data.frame()
  nonAvailableTiles <- predictionFrequencies %>% dplyr::filter(Freq == 0, Prediction != Reference)

  nonAvailableTiles$Domain <- NA
  nonAvailableTiles$Prediction <- as.character(nonAvailableTiles$Prediction)
  nonAvailableTiles$Reference <- as.character(nonAvailableTiles$Reference)

  for (i in seq(1:nrow(nonAvailableTiles))) {

    if (nonAvailableTiles$Prediction[i] %in% abbreviations[, classColumn]) {
      nonAvailableTiles$Prediction[i] <- abbreviations[abbreviations[,classColumn] == nonAvailableTiles$Prediction[i], "abbreviation"] %>% base::unique()
    }
    if (nonAvailableTiles$Reference[i] %in% abbreviations[, classColumn]) {
      nonAvailableTiles$Reference[i] <- abbreviations[abbreviations[, classColumn] == nonAvailableTiles$Reference[i], "abbreviation"]  %>% base::unique()
    }
  }
  return(nonAvailableTiles)

}
