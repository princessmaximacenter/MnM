integrateScores <- function(resultHigherClass) {
  resultHigherClass$finalHigherClass <- NA
  resultHigherClass$finalHigherClassScore <- 0
  for (i in seq(1:nrow(resultHigherClass))) {
    if (resultHigherClass$higherClassScore[i] >= resultHigherClass$higherClassScore2[i]) {
      resultHigherClass$finalHigherClass[i] <- resultHigherClass$predictHigherClass[i]
      resultHigherClass$finalHigherClassScore[i] <- resultHigherClass$higherClassScore[i]
    } else if (resultHigherClass$higherClassScore2[i] > resultHigherClass$higherClassScore[i]) {
      resultHigherClass$finalHigherClass[i] <- resultHigherClass$predictHigherClass2[i]
      resultHigherClass$finalHigherClassScore[i] <- resultHigherClass$higherClassScore2[i]
    }
  }
  return(resultHigherClass)
}
