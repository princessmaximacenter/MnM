#' Calcalate metrics to estimate probability score treshold
#'
#' This function calculates the sensitivity, specificity, precision and recall at
#' different probability score thresholds for sample classification.
#' A sample is classified as positive when it receives a classification (prediction probability > threshold).
#' On the contrary, a sample is classified as negative when it does not receive a classification
#' (prediction probability < threshold). A true positive is a classified sample for which the prediction is
#' correct, a true negative a non-classified sample that had a wrong prediction.
#'
#' @param predictionsMM  Dataframe showing the top 3 classifications for the tumor (sub)type,
#' together with their probability scores and the original pathology label (originalCall)
#'
#' @return Dataframe containing the number of true positives ($TP), false positives ($FP),
#' true negatives ($TN) and false negatives ($FN) at different probability score tresholds.
#'  Calculated from the metrics named above are the sensitivity ($sensitivity),
#'  specificity ($specificity), precision ($precision) and recall ($recall).
#' @export
#'
#'

estimateProbabilityThreshold <- function(predictionsMM) {

  for (cutoff in seq(from = 0, to = 1.01, by = 0.01)) {
    predictionsMMFiltered <- predictionsMM %>% filter(probability1 > cutoff)

    TP <- predictionsMMFiltered %>% filter(originalCall == predict) %>% nrow()
    FP <- predictionsMMFiltered %>% filter(originalCall != predict) %>% nrow()
    TN <- predictionsMM %>% filter(probability1 <= cutoff,
                                   originalCall != predict) %>% nrow()
    FN <- predictionsMM %>% filter(probability1 <= cutoff,
                                   originalCall == predict) %>% nrow()

    #recall <- nrow(predictionsMMFiltered) / nrow(predictionsMM)
    binaryScoreDF <- data.frame(cutoff = cutoff,
                                TP = TP,
                                FP = FP,
                                TN = TN,
                                FN = FN)
    if (cutoff == 0) {
      totalbinaryScoreDF <- binaryScoreDF
    } else {
      totalbinaryScoreDF <- rbind(totalbinaryScoreDF,
                                  binaryScoreDF)
    }
  }
  totalbinaryScoreDF$sensitivity <- totalbinaryScoreDF$TP/ (totalbinaryScoreDF$TP + totalbinaryScoreDF$FN)

  totalbinaryScoreDF$specificity <- totalbinaryScoreDF$TN/ (totalbinaryScoreDF$TN + totalbinaryScoreDF$FP)

  totalbinaryScoreDF$precision <- totalbinaryScoreDF$TP / (totalbinaryScoreDF$TP + totalbinaryScoreDF$FP)

  totalbinaryScoreDF$recall <- (totalbinaryScoreDF$TP + totalbinaryScoreDF$FP) /
    (totalbinaryScoreDF$TP + totalbinaryScoreDF$FP + totalbinaryScoreDF$TN + totalbinaryScoreDF$FN)

  totalbinaryScoreDF$precision[is.na(totalbinaryScoreDF$precision)] <- 1
  return(totalbinaryScoreDF)
}
