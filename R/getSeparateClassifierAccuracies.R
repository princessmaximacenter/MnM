#' Get separate results from minority and majority classifier
#'
#' @param minority R-object that contains the results from the Minority classifier
#' @param majority R-object that contains the results from the Majority classifier
#' @param crossValidation  Specify whether the results are from the cross-validation setup or not.
#' This is important, as for the cross-validation setup there is a ground truth ($originalCall), while for new predictions there is not.
#' @param classColumn Column name within metadata-file that contains the cancer subtype-labels.
#' @param higherClassColumn Column name within metadata-file that contains the cancer type labels.
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' @param nModels How many models were used to obtain a final prediction?
#' @param probabilityThreshold What is the probability score threshold that needs to be used to filter confident sample classifications with?
#' @return Dataframe containing the total amount of samples ($nSamples), how many samples pass the the probability score threshold ($nSamplesFiltered),
#' what is the fraction non-classified samples ($notClassified),
#' what is the fraction of samples that pass the threshold ($fraction), what is the fraction of correct samples within the highest scoring classification ($fractionCorrect),
#' the top 2 scoring classification ($fractionCorrect2) and top 3 scoring classifications ($fractionCorrect3), what is the fraction of incorrect classifications ($fractionIncorrect),
#' how many of the confident classifications are correct ($fractionCorrectFiltered), and certain score metrics ($Precision, $F1 and $Recall),
#' all split up based on their population frequency within the reference cohohrt ($nCases) and the classifier used to obtain the metrics.
#' @export
#'
getSeparateClassifierAccuracies <- function(minority,
                                             majority,
                                             crossValidation,
                                             classColumn,
                                             higherClassColumn,
                                             subtype = F,
                                            probabilityThresholdMajority,
                                            probabilityThresholdMinority,
                                             nModels) {

  predictionsList <- integrateMM(minority = minority,
              majority = majority,
              metaDataRef = metaDataRef,
              nModels = nModels,
              subtype = subtype,
              classColumn = classColumn,
              higherClassColumn = higherClassColumn,
              crossValidation = crossValidation,
              integrate = F)

  predictionsMajority <- predictionsList$predictionsMajority
  predictionsMinority <- predictionsList$predictionsMinority

  if (subtype == F) {

    minorityAccuracies <- getAccuraciesPerTumorTypeSize(predictionsMM = predictionsMinority,
                                                        metaDataRef = minority$metaData,
                                                        classColumn = higherClassColumn,
                                                        probabilityThreshold = probabilityThresholdMinority)

    majorityAccuracies <- getAccuraciesPerTumorTypeSize(predictionsMM = predictionsMajority,
                                                        metaDataRef = majority$metaData,
                                                        classColumn = higherClassColumn,
                                                        probabilityThreshold = probabilityThresholdMajority)
  } else {
    minorityAccuracies <- getAccuraciesPerTumorTypeSize(predictionsMM = predictionsMinority,
                                                        metaDataRef = minority$metaData,
                                                        classColumn = classColumn,
                                                        probabilityThreshold = probabilityThresholdMinority)

    majorityAccuracies <- getAccuraciesPerTumorTypeSize(predictionsMM = predictionsMajority,
                                                        metaDataRef = majority$metaData,
                                                        classColumn = classColumn,
                                                        probabilityThreshold = probabilityThresholdMajority)

  }

  minorityAccuracies$classifier <- "Minority Classifier"

  majorityAccuracies$classifier <- "Majority Classifier"

  fractionsCorrectTotal <- rbind(minorityAccuracies, majorityAccuracies)

  return(fractionsCorrectTotal)
}
