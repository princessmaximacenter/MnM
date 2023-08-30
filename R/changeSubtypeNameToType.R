#' Function to find the associated tumor type for each tumor subtype prediction.
#'
#' @param probabilities probabilities from the different predictions from the minority or majority classifier.
#' @param linkClassAndHigherClass dataframe containing the link between the tumor subtypes and the tumor types.
#' The dataframe is automatically generated in previous functions from the metadata.
#'
#' @return
#' @export
#'
#' @examples
changeSubtypeNameToType <- function(probabilities, linkClassAndHigherClass) {

  for (i in seq(1:length(probabilities))) {

    allNames <- names(probabilities[[i]]) %>% unlist(., use.names = F)

    #firstTry <- unlist(firstTry, use.names = F)
    allNamesFinal <- c()

    for (j in seq(1:length(allNames))) {
      allNamesFinal[j] <- linkClassAndHigherClass[allNames[j] == linkClassAndHigherClass$Disease_sub_specification1, "Disease_sub_class"]
    }
    names(probabilities[[i]]) <- allNamesFinal
  }

  # Sum the predictions that lead to the same tumor type
  probabilitiesFinal <- probabilities
  for (i in seq(1:length(probabilities))) {
    probabilitiesFinal[[i]] <- tapply(probabilities[[i]], names(probabilities[[i]]), sum)
  }

  return(probabilitiesFinal)
}
