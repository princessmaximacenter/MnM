#' Integrate results Minority & Majority classifier
#'
#' Function to obtain final classifications for the M&M algorithm by combining the results of the Minority & Majority classifiers.
#'
#' @param minority R-object that contains the results from the Minority classifier
#' @param majority R-object that contains the results from the Majority classifier
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level (TRUE) or
#' tumor type level (FALSE)?
#' @param integrate Do you want to integrate the results of the Minority and Majority classifiers into M&M,
#' or do you want to obtain the results from the separate classifiers?
#' If so, use _subtype = T_. If you want to obtain tumor type classifications instead, use _subtype = F_.
#' @return If integrate = T: List containing the a dataframe with the final classifications (predictionsMMFinal),
#' and a list with the probability scores for all classification labels that were assigned to samples (MMProbabilityList).
#'
#' predictionsMMFinal contains the top 3 final classification labels ($predict{2,3}) with their accompanying probability scores
#' ($probability{1,2,3} and the original diagnosis label ($originalCall).
#'
#' If integrate = F: List containing a dataframe with the final classifications separately for the Minority Classifier ($predictionsMinority)
#' and Majority Classifier ($predictionsMajority). This dataframe has the same setup as predictionsMMFinal.
#'
#' @export
#'
integrateMM <- function(minority,
                        majority,
                        subtype,
                        integrate = T
                        ) {

  classColumn <- minority$metaDataRun$classColumn
  higherClassColumn <- minority$metaDataRun$higherClassColumn
  probabilitiesMinority <- obtainProbabilities(minority)
  probabilitiesMajority <- obtainProbabilities(majority)

  if (base::length(probabilitiesMinority) > base::length(probabilitiesMajority)) {
    probabilitiesMinority <- probabilitiesMinority[base::names(probabilitiesMajority)]
  } else if (length(probabilitiesMinority) < base::length(probabilitiesMajority)) {
    probabilitiesMajority <- probabilitiesMajority[base::names(probabilitiesMinority)]
  }
  if (subtype == F) {
    linkClassAndHigherClass <- minority$metaDataRef[ , c(classColumn, higherClassColumn)] %>% base::unique()

    probabilitiesMinority <- changeSubtypeNameToType(probabilitiesMinority,
                                                          linkClassAndHigherClass = linkClassAndHigherClass,
                                                     classColumn = classColumn,
                                                     higherClassColumn = higherClassColumn)
    probabilitiesMajority <- changeSubtypeNameToType(probabilitiesMajority,
                                                          linkClassAndHigherClass = linkClassAndHigherClass,
                                                     classColumn = classColumn,
                                                     higherClassColumn = higherClassColumn)
  }
  if (integrate == T) {
  MMProbabilityList <- getMMProbabilities(minorityProbability = probabilitiesMinority,
                                          majorityProbability = probabilitiesMajority)

  predictionsMM <- getTopClassifications(minority = minority,
                                          MMProbabilityList = MMProbabilityList,
                                          higherClassColumn = higherClassColumn,
                                          subtype = subtype)

  predictionsList <- list(predictionsMMFinal = predictionsMM,
                          MMProbabilityList = MMProbabilityList)
  } else {
    predictionsMinority <- getTopClassifications(minority = minority,
                                                  MMProbabilityList = probabilitiesMinority,
                                                  higherClassColumn = higherClassColumn,
                                                  subtype = subtype)

    predictionsMajority <- getTopClassifications(minority = majority,
                                                  MMProbabilityList = probabilitiesMajority,
                                                  higherClassColumn = higherClassColumn,
                                                  subtype = subtype)

    predictionsList <- list(predictionsMinority = predictionsMinority,
                            predictionsMajority = predictionsMajority)
  }
  return(predictionsList)
}
