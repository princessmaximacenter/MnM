#' Check out the F1 score metrics for whole dataset
#'
#' @param predictionsMM Final predictions with the top 3 final predictions ($predict{2,3})
#' with their accompanying probability scores ($probability{1,2,3} and the original diagnosis label ($originalCall).
#' @return Dataframe containing the macro-F1 score for the unfiltered (macroF1) and filtered dataset (macroF1FilteredWithNA, macroF1Filtered), median F1 score
#' for the unfiltered (medianF1) and filtered dataset (medianF1Filtered).
#'
getF1ScoreMetrics <- function(predictionsMM) {

  tumorConfusionMatrix <- confusionMatrix(factor(predictionsMM$predict,
                                                 levels = unique(predictionsMM$originalCall)),
                                          factor(predictionsMM$originalCall, levels = unique(predictionsMM$originalCall)),
                                          dnn = c("Prediction", "Reference"))

  confusionDF <- tumorConfusionMatrix$byClass %>% as.data.frame()

  macroF1 <- mean(confusionDF$F1, na.rm = T)
  medianF1 <- median(confusionDF$F1, na.rm = T)

  predictionsMMFiltered <- predictionsMM %>% filter(probability1 > 0.8)

  tumorConfusionMatrixFiltered <- confusionMatrix(factor(predictionsMMFiltered$predict,
                                                         levels = unique(c(predictionsMMFiltered$originalCall, predictionsMMFiltered$predict))),
                                                  factor(predictionsMMFiltered$originalCall, levels = unique(c(predictionsMMFiltered$originalCall,
                                                                                                               predictionsMMFiltered$predict))),
                                                  dnn = c("Prediction", "Reference"))

  confusionDFFiltered <- tumorConfusionMatrixFiltered$byClass %>% as.data.frame()

  macroF1FilteredWithNA <- mean(confusionDFFiltered$F1, na.rm = T)
  confusionDFFiltered[is.na(confusionDFFiltered$F1),"F1"] <- 0
  macroF1Filtered <- mean(confusionDFFiltered$F1, na.rm = T)
  medianF1Filtered <- median(confusionDFFiltered$F1, na.rm = T)
  F1ScoreDF <- data.frame(macroF1 = macroF1,
                          medianF1 = medianF1,
             macroF1FilteredWithNA = macroF1FilteredWithNA,
             macroF1Filtered = macroF1Filtered,
             medianF1Filtered = medianF1Filtered)

  return(F1ScoreDF)
}
