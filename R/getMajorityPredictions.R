#' Obtain M&M majority predictions
#'
#' Function to extract the top 3 tumor (sub)type predictions for samples, with their accompanying probability scores.
#' @param minority R-object that contains the results from the Minority classifier.
#' @param MMProbabilityList List of all the samples containing the probabilities for the different sample predictions from the integrated M&M classifier.
#' @param higherClassColumn Column name within metadata-file that contains the cancer type labels.
#' @param crossValidation Specify whether the results are from the cross-validation setup or not.
#' @param metaDataRef Metadata-file for the reference cohort.
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' If so, use _subtype = T_. If you want to obtain tumor type predictions instead, use _subtype = F_.
#'
#' @return Dataframe showing the top 3 predictions for the tumor (sub)type, together with their probability scores.
#' @export
#'
getMajorityPredictions <- function(minority,
                                   MMProbabilityList,
                                   higherClassColumn,
                                   crossValidation,
                                   metaDataRef,
                                   subtype) {

  if (crossValidation == T) {
  predictions <- minority$classifications[, c("predict", "originalCall")]
  if(subtype == F) {
    predictions$originalCall <- metaDataRef[rownames(predictions), higherClassColumn]
  }
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


  return(predictions)
}
