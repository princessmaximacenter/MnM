#' Link tumor subtype to tumor type
#'
#' Function to find the associated tumor type for each tumor subtype prediction.
#'
#' @param probabilities List with the probabilities from the different tumor subtype predictions from the minority or majority classifier.
#' @param linkClassAndHigherClass dataframe containing the link between the tumor subtypes and the tumor types.
#' The dataframe is automatically generated in previous functions from the metadata.
#'
#' @return  List of all the samples containing the probabilities for the different tumor type predictions from the minority or majority classifier.
#' @export
#'
changeSubtypeNameToType <- function(probabilities,
                                    linkClassAndHigherClass,
                                    classColumn,
                                    higherClassColumn
                                    ) {

  for (i in seq(1:length(probabilities))) {

    allNames <- names(probabilities[[i]]) %>% unlist(., use.names = F)

    #firstTry <- unlist(firstTry, use.names = F)
    allNamesFinal <- c()

    for (j in seq(1:length(allNames))) {
      allNamesFinal[j] <- linkClassAndHigherClass[allNames[j] == linkClassAndHigherClass[,classColumn], higherClassColumn]
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
