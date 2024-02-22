top3Predictions <- function(predictionsMM) {

  top2 <- predictionsMM %>% filter(originalCall == predict | originalCall == predict2) %>% nrow()
  top3 <- predictionsMM %>% filter(originalCall == predict | originalCall == predict2 |
                                     originalCall == predict3) %>% nrow()

  topOnes <- c("predict","predict2","predict3")

  totalPredictions <- nrow(predictionsMM)

  for (i in seq(1:length(topOnes))) {

    top <-  predictionsMM %>% filter(originalCall == !!sym(topOnes[i])) %>% nrow()

    if ( i == 1) {
      topper <- top
    } else {
      topper <- topper + top
    }

    topScores <- data.frame(topN = topOnes[i],
                            howManyClassifiedCorrectly = topper,
                            percentageClassifiedCorrectly = round(topper / totalPredictions,
                                                                  digits = 2))

    if (i == 1) {
      topScoresDF <- topScores
    } else {
      topScoresDF <- rbind(topScoresDF,
                           topScores)
    }

  }
  return(topScoresDF)

}

