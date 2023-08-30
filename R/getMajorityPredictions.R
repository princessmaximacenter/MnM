#' Title
#'
#' @param minority
#' @param MMProbabilityList
#' @param higherClassColumn
#' @param crossValidation
#'
#' @return
#' @export
#'
#' @examples
getMajorityPredictions <- function(minority,
                                   MMProbabilityList,
                                   higherClassColumn,
                                   crossValidation) {

  if (crossValidation == T) {
  predictions <- minority$classifications[, c("predict", "originalCall")]
  } else {
    predictions <- minority$classifications[, "predict", drop = F]
  }

  predictions$probability1 <- NA
  predictions$probability2 <- NA
  predictions$probability3 <- NA
  predictions$predict2 <- NA
  predictions$predict3 <- NA

  for (i in seq(1:length(MMProbabilityList))) {
    allProbs <- MMProbabilityList[[i]]
    numbersProbs <- as.vector(allProbs)
    names(numbersProbs) <- names(allProbs)
    numbersProbs <- sort(numbersProbs, decreasing = T)
    highestProbability <- numbersProbs[1]
    secondProbability <- numbersProbs[2]
    thirdProbability <- numbersProbs[3]

    predictions$probability1[i] <- highestProbability
    predictions$probability2[i] <- secondProbability
    predictions$probability3[i] <- thirdProbability

    predictions$predict[i] <- names(highestProbability)
    predictions$predict2[i] <- names(secondProbability)
    predictions$predict3[i] <- names(thirdProbability)
  }

  predictions$originalCall <- metaDataRef[rownames(predictions), higherClassColumn]

    notMalignant <- c("Not malignant bone marrow",
                      "Not malignant blood",
                      "Bone marrow failure" )

    gsubName <- "Not malignant (Hemato)"

  for (i in seq(1:length(notMalignant))) {

    predictions$originalCall <- gsub(notMalignant[i], gsubName, predictions$originalCall)

  }



  return(predictions)
}
