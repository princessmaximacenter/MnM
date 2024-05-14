#' Obtain probabilities for classifications
#'
#' This function calculates the probability for each classification,
#' by dividing the number of models for a specific tumor (sub)type label by the total number of models.
#' In this way, the probability is equal to the fraction of classifier calls for a certain label.
#'
#' @param classifierResults R-object containing the results from the Minority or Majority classifier.
#' The object should contain a list with the probabilities for each classification label

#' @return List of all the samples containing the probabilities for the different sample classifications from the generated models.
#'
obtainProbabilities <- function(classifierResults) {

  depth <- function(this) ifelse(is.list(this), 1L + base::max(base::sapply(this, depth)), 0L)
  probabilityList <- classifierResults$probabilityList

  if (depth(probabilityList) == 2) {
    for (i in base::seq(1:base::length(probabilityList))) {
      nModels <- base::sum(probabilityList[[1]][[1]])
      if (i == 1) {

        probabilities <- base::lapply(probabilityList[[i]], function(x) x / nModels)
      } else {
        probabilities <- c(probabilities, base::lapply(probabilityList[[i]], function(x) x / nModels))
      }
    }
  } else {
    nModels <- base::sum(probabilityList[[1]])

    probabilities <- base::lapply(classifierResults$probabilityList, function(x) x / nModels)
  }

  return(probabilities)
}
