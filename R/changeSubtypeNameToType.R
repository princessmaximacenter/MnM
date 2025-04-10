#' Link tumor subtype to tumor type
#'
#' Function to find the associated tumor type for each tumor subtype prediction.
#'
#' @param probabilities List with the probabilities from the different tumor subtype predictions from the Minority or Majority classifier.
#' @param linkClassAndHigherClass dataframe containing the link between the tumor subtypes and the tumor types.
#' The dataframe is automatically generated in previous functions from the metadata.
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Name of column in the metadata file that contains the tumor type labels.
#'
#' @return  List of all the samples containing the probabilities for the different _tumor type_ predictions from the Minority or Majority classifier.
#'
changeSubtypeNameToType <- function(probabilities,
                                    linkClassAndHigherClass,
                                    classColumn,
                                    higherClassColumn
                                    ) {

  for (i in base::seq(1:base::length(probabilities))) {

    allNames <- base::names(probabilities[[i]]) %>% base::unlist(., use.names = F)

    #firstTry <- unlist(firstTry, use.names = F)
    allNamesFinal <- c()

    for (j in base::seq(1:base::length(allNames))) {
      allNamesFinal[j] <- linkClassAndHigherClass[allNames[j] == linkClassAndHigherClass[,classColumn], higherClassColumn]
    }
    base::names(probabilities[[i]]) <- allNamesFinal
  }

  # Sum the predictions that lead to the same tumor type
  probabilitiesFinal <- probabilities
  for (i in base::seq(1:base::length(probabilities))) {
    probabilitiesFinal[[i]] <- base::tapply(probabilities[[i]], base::names(probabilities[[i]]), base::sum)
  }

  return(probabilitiesFinal)
}
