obtainProbabilityDF <- function(classifierResults, probabilities) {

  predictions <- classifierResults$classifications

  predictions$probability1 <- NA
  predictions$probability2 <- NA
  predictions$probability3 <- NA
  predictions$predict2 <- NA
  predictions$predict3 <- NA

  for (i in seq(1:length(probabilities))) {
    allProbs <- probabilities[[i]]
    numbersProbs <- as.vector(allProbs)
    names(numbersProbs) <- dimnames(allProbs)[[1]]
    numbersProbs <- sort(numbersProbs, decreasing = T)
    highestProbability <- numbersProbs[1]
    secondProbability <- numbersProbs[2]
    thirdProbability <- numbersProbs[3]

    predictions$probability1[i] <- highestProbability
    predictions$probability2[i] <- secondProbability
    predictions$probability3[i] <- thirdProbability

    predictions$predict2[i] <- names(secondProbability)
    predictions$predict3[i] <- names(thirdProbability)
  }

  return(predictions)

}
