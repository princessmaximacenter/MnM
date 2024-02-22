#' Calculate accuracy for score blocks combined
#'
#' This function calculates for all samples with a higher probability score than a certain treshold what is the accuracy.
#' Furthermore, the total number of samples that has a probability score higher than the threshold is determined.
#'
#' @param totalAccDF Dataframe showing what is the corresponding score block ($Filter),
#' how many samples fall within the probability score range ($total),
#' how many of these samples receive a correct prediction ($totalCorrect),
#' how many an incorrect prediction ($totalIncorrect),
#' and what is the accuracy within the score block ($accuracy).
#'
#' @return The same dataframe as before, but now containing new columns $cumulativeAccuracy
#' (combined accuracy also with the higher score blocks),
#' and $cumulativeTotal (how many samples are classified in the current and higher score blocks).
#'
getCumulativeAccuracy <- function(totalAccDF, totalSamples) {
  totalAccDF$cumulativeAccuracy  <- 0
  totalAccDF$cumulativeTotal <- 0

  for (i in seq(1:nrow(totalAccDF))) {
    if (i == 1) {
      total <- totalAccDF$total[i]
      totalCorrect <- totalAccDF$totalCorrect[i]
      totalAccDF$cumulativeAccuracy[i] <- totalAccDF$accuracy[i]
      totalAccDF$cumulativeTotal[i] <- total
      totalAccDF$percentClassified[i] <- round(total / totalSamples, digits = 2)
    } else {
      total <- total + totalAccDF$total[i]
      totalCorrect <- totalCorrect + totalAccDF$totalCorrect[i]
      totalAccDF$cumulativeAccuracy[i] <- round(totalCorrect/total * 100, digits = 1)
      totalAccDF$cumulativeTotal[i] <- total
      totalAccDF$percentClassified[i] <- round(total / totalSamples, digits = 3)
    }

  }
  return(totalAccDF)
}
