#' Integrate results Minority & Majority classifier
#'
#' Function to obtain final predictions for the M&M algorithm by combining the results of the Minority & Majority classifiers.
#'
#' @param minority R-object that contains the results from the Minority classifier
#' @param majority R-object that contains the results from the Majority classifier
#' @param metaDataRef Metadata-file for the reference cohort.
#' @param nModels How many models were used to obtain a final prediction?
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' If so, use _subtype = T_. If you want to obtain tumor type predictions instead, use _subtype = F_.
#' @param classColumn Column name within metadata-file that contains the cancer subtype-labels.
#' @param higherClassColumn Column name within metadata-file that contains the cancer type labels.
#' @param crossValidation Specify whether the results are from the cross-validation setup or not.
#' This is important, as for the cross-validation setup there is a ground truth ($originalCall), while for new predictions there is not.
#'
#' @return Dataframe showing the top 3 predictions for the tumor (sub)type, together with their probability scores.
#' @export
#'
integrateMM <- function(minority,
                        majority,
                        metaDataRef,
                        nModels,
                        subtype,
                        classColumn,
                        higherClassColumn,
                        crossValidation = T) {

  probabilitiesMinority <- obtainProbabilities(minority,
                                               crossValidation = crossValidation,
                                               nModels = nModels
                                               )
  probabilitiesMajority <- obtainProbabilities(majority,
                                               crossValidation = crossValidation,
                                               nModels = nModels)


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
  MMProbabilityList <- getMMProbabilities(minorityProbability = probabilitiesMinority,
                                          majorityProbability = probabilitiesMajority)

  predictionsMM <- getMajorityPredictions(minority = minority,
                                          MMProbabilityList = MMProbabilityList,
                                          higherClassColumn = higherClassColumn,
                                          crossValidation = crossValidation,
                                          metaDataRef = metaDataRef,
                                          subtype = subtype)

  return(predictionsMM)
}
