#' Obtain M&M's top classifications
#'
#' Function to extract the top 3 tumor (sub)type classifications for samples, with their accompanying probability scores.
#' @param minority R-object that contains the results from the Minority classifier.
#' @param MMProbabilityList List of all the samples containing the probabilities for the different sample predictions from the integrated M&M classifier.
#' @param higherClassColumn Column name within metadata-file that contains the tumor type labels.
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' If so, use _subtype = T_. If you want to obtain tumor type predictions instead, use _subtype = F_.
#'
#' @return Dataframe showing the top 3 classifications for the tumor (sub)type ($predict{2,3}), together with their probability scores ($probability{1,2,3}.
#' @export
#'
getTopClassifications <- function(minority,
                                   MMProbabilityList,
                                   higherClassColumn,
                                   subtype) {


  if ("originalCall" %in% colnames(minority$classifications)) {
  predictions <- minority$classifications[base::names(MMProbabilityList), c("predict", "originalCall")]
  if(subtype == F) {
    predictions$originalCall <- minority$metaDataRef[rownames(predictions), higherClassColumn]
  }
  } else {
    predictions <- minority$classifications[, "predict", drop = F]
  }

  predictions$probability1 <- NA
  predictions$predict2 <- NA
  predictions$probability2 <- NA
  predictions$predict3 <- NA
  predictions$probability3 <- NA

  for (i in seq(1:base::length(MMProbabilityList))) {
    allProbs <- MMProbabilityList[[i]]
    numbersProbs <- as.vector(allProbs)
    names(numbersProbs) <- base::names(allProbs)
    numbersProbs <- base::sort(numbersProbs, decreasing = T)
    highestProbability <- numbersProbs[1]
    secondProbability <- numbersProbs[2]
    thirdProbability <- numbersProbs[3]

    predictions$probability1[i] <- highestProbability
    predictions$probability2[i] <- secondProbability
    predictions$probability3[i] <- thirdProbability

    predictions$predict[i] <- base::names(highestProbability)
    predictions$predict2[i] <- base::names(secondProbability)
    predictions$predict3[i] <- base::names(thirdProbability)
  }


  return(predictions)
}
