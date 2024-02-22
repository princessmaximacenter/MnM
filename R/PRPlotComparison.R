PRPlotComparison <- function(predictionsMM, otherClassifierName,
                             otherClassifierResults,
                             originalCallColumnOtherClassifier,
                             otherClassifierPredictionColumn,
                             scoreName) {
otherClassifierResultsMnM <- predictionsMM %>% filter(rownames(.) %in% rownames(otherClassifierResults))
ROCcurveMnM <- estimateProbabilityThreshold(otherClassifierResultsMnM)

for (cutoff in seq(from = 0, to = 1.01, by = 0.01)) {
  predictionsFiltered <- otherClassifierResults %>% filter(!!sym(scoreName) >= cutoff)

  TP <- predictionsFiltered %>% filter(!!sym(originalCallColumnOtherClassifier) == !!sym(otherClassifierPredictionColumn)) %>% nrow()
  FP <- predictionsFiltered %>% filter(!!sym(originalCallColumnOtherClassifier) != !!sym(otherClassifierPredictionColumn)) %>% nrow()
  TN <- otherClassifierResults %>% filter(!!sym(scoreName) < cutoff,
                                !!sym(originalCallColumnOtherClassifier) != !!sym(otherClassifierPredictionColumn)) %>% nrow()
  FN <- otherClassifierResults %>% filter(!!sym(scoreName) < cutoff,
                                !!sym(originalCallColumnOtherClassifier) == !!sym(otherClassifierPredictionColumn)) %>% nrow()
  #recall <- nrow(predictionsFiltered) / nrow(otherClassifierResults)
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

totalbinaryScoreDF$type <- otherClassifierName
ROCcurveMnM$type <- "M&M"

totalDF <- rbind(totalbinaryScoreDF,
                 ROCcurveMnM)

colnames(totalDF)[6:10] <- c("Sensitivity", "Specificity",
                             "Precision", "Recall", "type")

#newlineMnM <- data.frame(cutoff = NA,
#                          TP = NA,
#                          FP = NA,
#                          TN = NA,
#                          FN = NA,
#                          Sensitivity = NA,
#                          Specificity = NA,
#                          Precision = 0,
#                          Recall = 1,
#                          type = "M&M")
#
# newlineOther <- data.frame(cutoff = NA,
#                          TP = NA,
#                          FP = NA,
#                          TN = NA,
#                          FN = NA,
#                          Sensitivity = NA,
#                          Specificity = NA,
#                          Precision = 0,
#                          Recall = 1,
#                          type = otherClassifierName)
#totalDF <- rbind(newlineOther,
#      newlineMnM,
#      totalDF)

return(totalDF)
}
