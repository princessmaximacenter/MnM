#' Determine link accuracy and probability score
#'
#' Get the accuracies of predictions associated with specific probability scores.
#' Score blocks will have a range of 0.1, starting from 1 and incrementally decreasing to 0, or higher if specified with the _sequencesMax_ argument.
#'
#' @param predictionsMM Dataframe with the predictions from the M&M classification process.
#' @param sequencesMax Do you want to stop at a specified probability score threshold? If so, specify here.
#'
#' @return Dataframe showing what is the corresponding score block ($Filter), how many samples fall within the probability score range ($total),
#' how many of these samples receive a correct prediction ($totalCorrect), how many an incorrect prediction ($totalIncorrect),
#' what is the accuracy within the score block ($accuracy), what is the combined accuracy also with the higher score blocks ($cumulativeAccuracy),
#' and how many samples are classified in the current and higher score blocks ($cumulativeTotal).
#' @export
#'
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
  totalAccDF <- getCumulativeAccuracy(totalAccDF, totalSamples = nrow(predictionsMM))

  return(totalAccDF)
}
