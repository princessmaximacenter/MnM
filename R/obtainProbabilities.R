#' Obtain probabilities for classifier predictions
#'
#' This function calculates the probability for each prediction, by dividing the number of models predicting a tumor (sub)type by the total number of models ( _nModels_).
#'
#' @param classifierResults R-object containing the results from the minority or majority classifier.
#' The object should contain a list with the probabilities for each prediction.
#' @param crossValidation A boolean factor specifying whether the R-object was generated within a cross-validation setup.
#' @param nModels How many models were used to obtain a final prediction?
#'
#' @return List of all the samples containing the probabilities for the different sample predictions from the _nModels_ generated models.
#'
obtainProbabilities <- function(classifierResults,
                                   crossValidation,
                                   nModels) {

  if (crossValidation == T) {
    probabilityList <- classifierResults$probabilityList

    for (i in seq(1:length(probabilityList))) {
      if (i == 1) {
        probabilities <- lapply(probabilityList[[i]], function(x) x / nModels)
      } else {
        probabilities <- c(probabilities, lapply(probabilityList[[i]], function(x) x / nModels))
      }
    }
  } else {
    probabilities <- lapply(classifierResults$probability, function(x) x / nModels)
  }

  return(probabilities)
}
