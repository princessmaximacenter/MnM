getAccuraciesPerScoreBlock <- function(predictionsMM,
                                               sequencesMax = NA) {

  if (is.na(sequencesMax)) {
    sequenceHigher <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0)
    sequenceLower <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
  } else {
    sequenceHigher <- seq(0.9, sequencesMax, by = -0.1)
    sequenceLower <- seq(1, sequencesMax + 0.1, by = -0.1)
  }

  for (i in seq(1:length(sequenceHigher))) {

    totalCorrect <- predictionsMM %>% filter(probability1 > sequenceHigher[i] & probability1 <= sequenceLower[i],
                                                  predict == originalCall) %>%
      nrow()

    totalIncorrect <- predictionsMM %>% filter(probability1 > sequenceHigher[i] & probability1 <= sequenceLower[i],
                                                    predict != originalCall) %>%
      nrow()

    accuracyDF <- data.frame(Type = "M&M", Filter = sequenceHigher[i],
                             total = totalCorrect + totalIncorrect,
                             totalCorrect = totalCorrect,
                             totalIncorrect = totalIncorrect,
                             accuracy = round(totalCorrect/(totalCorrect + totalIncorrect) * 100, digits = 1))

    if (i == 1) {
      totalAccDF <- accuracyDF

    } else {
      totalAccDF <- rbind(totalAccDF, accuracyDF )
    }
  }
  totalAccDF <- getCumulativeAccuracy(totalAccDF)

  return(totalAccDF)
}
