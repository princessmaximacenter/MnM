#' Integrate Minority & Majority probability scores
#'
#' This function calculates the probability for each prediction,
#' summing up the probabilities for the Minority and Majority classifiers and dividing the probabilities either
#' by 2 (when both classifiers made the prediction) or 10 (when only one classifier made the prediction).
#'
#' @param majorityProbability List of all samples containing the probabilities for the different sample predictions from the Majority classifier.
#' @param minorityProbability List of all samples containing the probabilities for the different sample predictions from the Minority classifier.
#'
#' @return List of all the samples containing the probabilities for the different sample predictions from the integrated M&M classifier.
#' @export
#'
getMMProbabilities <- function(majorityProbability,
                               minorityProbability) {

  MMProbabilityList <- majorityProbability
  for (i in seq(1:length(minorityProbability))) {

    combineThem <- c(minorityProbability[[i]], majorityProbability[[names(minorityProbability[i])]])
    combineThem <- combineThem[order(names(combineThem), decreasing = T)]

    combinedThings <- c()
    tumors <- c()

    for (item in seq(from = 1,   to = (length(combineThem) - 1))) {
      nextItem <- item + 1
      if ((names(combineThem)[item] == names(combineThem)[nextItem]) &
          (names(combineThem)[item] %notin% tumors)) {
        combinedThings <- c(combinedThings, ((combineThem[item] + combineThem[nextItem]) / 2))
        tumors <- c(tumors, names(combineThem)[item])
      } else {
        if (names(combineThem)[item] %notin% tumors) {
          combinedThings <- c(combinedThings, combineThem[item] * 0.1 ) # Je wilt dat deze score laag is omdat dit maar door 1 van de 2 wordt gezegd
          tumors <- c(tumors, names(combineThem)[item])
        } else if ((names(combineThem)[nextItem] %notin% tumors ) &
                   nextItem == length(combineThem)) {
          combinedThings <- c(combinedThings, combineThem[nextItem] * 0.1)
        }
      }
    }
    MMProbabilityList[[names(minorityProbability[i])]] <- combinedThings
  }

  return(MMProbabilityList)

}
