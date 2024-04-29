#' Integrate Minority & Majority probability scores
#'
#' This function calculates the probability for each prediction,
#' summing up the probabilities for the Minority and Majority classifiers and dividing the probabilities either
#' by 2 (when both classifiers made the prediction) or 10 (when only one classifier made the prediction).
#'
#' @param majorityProbability List of all samples containing the probabilities for the different sample predictions from the Majority classifier.
#' @param minorityProbability List of all samples containing the probabilities for the different sample predictions from the Minority classifier.
#'
#' @return List of all the samples containing the probabilities for the different sample classifications from the integrated M&M classifier.
#'
getMMProbabilities <- function(majorityProbability,
                               minorityProbability) {

  `%notin%` <- Negate(`%in%`)
  MMProbabilityList <- majorityProbability

  if (base::length(minorityProbability) != base::length(majorityProbability)) {
    stop("There are different numbers of classifications within the Minority and Majority classifier results. Please make sure the classifiers have been run on the same samples.")
  }
  for (i in seq(1:base::length(minorityProbability))) {

    combineThem <- c(minorityProbability[[i]], majorityProbability[[base::names(minorityProbability[i])]])
    combineThem <- combineThem[base::order(base::names(combineThem), decreasing = T)]

    combinedThings <- c()
    tumors <- c()

    for (item in seq(from = 1, to = (base::length(combineThem) - 1))) {
      nextItem <- item + 1
      if ((base::names(combineThem)[item] == base::names(combineThem)[nextItem]) &
          (base::names(combineThem)[item] %notin% tumors)) {
        combinedThings <- c(combinedThings, ((combineThem[item] + combineThem[nextItem]) / 2))
        tumors <- c(tumors, base::names(combineThem)[item])
      } else {
        if (names(combineThem)[item] %notin% tumors) {
          combinedThings <- c(combinedThings, combineThem[item] / 10 )
        }
        if (nextItem == base::length(combineThem) &
            base::names(combineThem)[nextItem] %notin% tumors) {
          combinedThings <- c(combinedThings, combineThem[nextItem] /10 )
          tumors <- c(tumors, base::names(combineThem)[nextItem])
        }
      }
    }
    MMProbabilityList[[base::names(minorityProbability[i])]] <- combinedThings
  }

  return(MMProbabilityList)

}
