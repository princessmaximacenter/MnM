estimateProbabilityThreshold <- function(predictionsMM) {

  for (cutoff in seq(from = 0.01, to = 1.00, by = 0.01)) {
    predictionsMMFiltered <- predictionsMM %>% filter(probability1 > cutoff)

    TP <- predictionsMMFiltered %>% filter(originalCall == predict) %>% nrow()
    FP <- predictionsMMFiltered %>% filter(originalCall != predict) %>% nrow()
    TN <- predictionsMM %>% filter(probability1 <= cutoff,
                                   originalCall != predict) %>% nrow()
    FN <- predictionsMM %>% filter(probability1 <= cutoff,
                                   originalCall == predict) %>% nrow()

    binaryScoreDF <- data.frame(cutoff = cutoff,
                                TP = TP,
                                FP = FP,
                                TN = TN,
                                FN = FN)
    if (cutoff == 0.01) {
      totalbinaryScoreDF <- binaryScoreDF
    } else {
      totalbinaryScoreDF <- rbind(totalbinaryScoreDF,
                                  binaryScoreDF)
    }
  }
  totalbinaryScoreDF$sensitivity <- totalbinaryScoreDF$TP/ (totalbinaryScoreDF$TP + totalbinaryScoreDF$FN)

  totalbinaryScoreDF$specificity <- totalbinaryScoreDF$TN/ (totalbinaryScoreDF$TN + totalbinaryScoreDF$FP)

  totalbinaryScoreDF$recall <- totalbinaryScoreDF$sensitivity

  totalbinaryScoreDF$precision <- totalbinaryScoreDF$TP / (totalbinaryScoreDF$TP + totalbinaryScoreDF$FP)


  totalbinaryScoreDFGgplot <- totalbinaryScoreDF %>% pivot_longer(cols = sensitivity:precision,
                                                                  values_to = "counts")

  return(totalbinaryScoreDFGgplot)
}
