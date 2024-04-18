#' Obtain probabilities for classifier predictions
#'
#' This function calculates the probability for each prediction, by dividing the number of models predicting a tumor (sub)type by the total number of models ( _nModels_).
#'
#' @param classifierResults R-object containing the results from the minority or majority classifier.
#' The object should contain a list with the probabilities for each prediction.

#' @return List of all the samples containing the probabilities for the different sample predictions from the _nModels_ generated models.
#'
obtainProbabilities <- function(classifierResults) {

  depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, depth)), 0L)
  probabilityList <- classifierResults$probabilityList

  if (depth(probabilityList) == 2) {
    for (i in seq(1:length(probabilityList))) {
      nModels <- sum(probabilityList[[1]][[1]])
      if (i == 1) {

        probabilities <- lapply(probabilityList[[i]], function(x) x / nModels)
      } else {
        probabilities <- c(probabilities, lapply(probabilityList[[i]], function(x) x / nModels))
      }
    }
  } else {
    nModels <- sum(probabilityList[[1]])
    print(paste0("Converting label counts to probability scores for the ", nModels, " models."))

    probabilities <- lapply(classifierResults$probabilityList, function(x) x / nModels)
  }

  return(probabilities)
}
