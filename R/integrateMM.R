#' Integrate results Minority & Majority classifier
#'
#' Function to obtain final predictions for the M&M algorithm by combining the results of the Minority & Majority classifiers.
#'
#' @param minority R-object that contains the results from the Minority classifier
#' @param majority R-object that contains the results from the Majority classifier
#' @param metaDataRef Metadata file containing the links between the patients and
#' the tumor (sub)type diagnosis within the reference cohort.
#' @param nModels How many models were used to obtain a final prediction?
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' If so, use _subtype = T_. If you want to obtain tumor type predictions instead, use _subtype = F_.
#' @param classColumn Column name within metadata-file that contains the cancer subtype-labels.
#' @param higherClassColumn Column name within metadata-file that contains the cancer type labels.
#' @param crossValidation Specify whether the results are from the cross-validation setup or not.
#' This is important, as for the cross-validation setup there is a ground truth ($originalCall), while for new predictions there is not.
#' @param substituteNames Do you want to substitute some names?
#' @param substituteBy What do you want to substitute it by?
#' @return If integrate = T: List containing the a dataframe with the final predictions (predictionsMMFinal),
#' and a list with the probability scores for all classification labels that were assigned to samples (MMProbabilityList).
#'
#' predictionsMMFinal contains the top 3 final classification labels ($predict{2,3}) with their accompanying probability scores
#' ($probability{1,2,3} and the original diagnosis label ($originalCall).
#'
#' If integrate = F: List containing a dataframe with the final predictions separately for the Minority Classifier ($predictionsMinority)
#' and Majority Classifier ($predictionsMajority)
#'
#' @export
#'
integrateMM <- function(minority,
                        majority,
                        metaDataRef,
                        nModels,
                        subtype,
                        classColumn,
                        higherClassColumn,
                        crossValidation = T,
                        integrate = T


                        ) {

  probabilitiesMinority <- obtainProbabilities(minority,
                                               crossValidation = crossValidation,
                                               nModels = nModels
                                               )
  probabilitiesMajority <- obtainProbabilities(majority,
                                               crossValidation = crossValidation,
                                               nModels = nModels)

  if (length(probabilitiesMinority) > length(probabilitiesMajority)) {
    probabilitiesMinority <- probabilitiesMinority[names(probabilitiesMajority)]
  } else if (length(probabilitiesMinority) < length(probabilitiesMajority)) {
    probabilitiesMajority <- probabilitiesMajority[names(probabilitiesMinority)]
  }
  if (subtype == F) {
    linkClassAndHigherClass <- metaDataRef[ , c(classColumn, higherClassColumn)] %>% unique

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

  predictionsMM <- getMajorityPredictions(minority = minority,
                                          MMProbabilityList = MMProbabilityList,
                                          higherClassColumn = higherClassColumn,
                                          crossValidation = crossValidation,
                                          metaDataRef = metaDataRef,
                                          subtype = subtype)

  predictionsList <- list(predictionsMMFinal = predictionsMM,
                          MMProbabilityList = MMProbabilityList)
  } else {
    predictionsMinority <- getMajorityPredictions(minority = minority,
                                                  MMProbabilityList = probabilitiesMinority,
                                                  higherClassColumn = higherClassColumn,
                                                  crossValidation = crossValidation,
                                                  metaDataRef = metaDataRef,
                                                  subtype = subtype)

    predictionsMajority <- getMajorityPredictions(minority = majority,
                                                  MMProbabilityList = probabilitiesMajority,
                                                  higherClassColumn = higherClassColumn,
                                                  crossValidation = crossValidation,
                                                  metaDataRef = metaDataRef,
                                                  subtype = subtype)

    predictionsList <- list(predictionsMinority = predictionsMinority,
                            predictionsMajority = predictionsMajority)
  }
  return(predictionsList)
}
