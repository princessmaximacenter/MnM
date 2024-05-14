#' Get separate results from Minority and Majority classifier
#'
#' @param minority R-object that contains the results from the Minority classifier
#' @param majority R-object that contains the results from the Majority classifier
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level (subtype = TRUE)?
#' @param probabilityThresholdMajority What is the probability score threshold that needs to be used to
#' filter confident sample classifications with within the Majority classifier?
#' @param probabilityThresholdMinority What is the probability score threshold
#' that needs to be used to filter confident sample classifications with within the Minority classifier?
#'
#' @return Dataframe containing the total amount of samples ($nSamples), how many samples pass the the probability score threshold ($nSamplesFiltered),
#' what is the fraction non-classified samples ($notClassified),
#' what is the fraction of samples that pass the threshold ($fraction), what is the fraction of correct samples within the highest scoring classification ($fractionCorrect),
#' the top 2 scoring classification ($fractionCorrect2) and top 3 scoring classifications ($fractionCorrect3), what is the fraction of incorrect classifications ($fractionIncorrect),
#' how many of the confident classifications are correct ($fractionCorrectFiltered), and certain score metrics ($Precision, $F1 and $Recall),
#' all split up based on their population frequency within the reference cohort ($nCases) and the classifier used to obtain the metrics
#' ($classifier = Minority / Majority).
#'
getSeparateClassifierAccuracies <- function(minority,
                                             majority,
                                             subtype = F,
                                            probabilityThresholdMajority = 0.9,
                                            probabilityThresholdMinority = 0.73
                                             ) {

  predictionsList <- integrateMM(minority = minority,
              majority = majority,
              subtype = subtype,
              integrate = F)

  predictionsMajority <- predictionsList$predictionsMajority
  predictionsMinority <- predictionsList$predictionsMinority
  classColumn <- minority$metaDataRun$classColumn
  higherClassColumn <- minority$metaDataRun$higherClassColumn
  if (subtype == F) {

    minorityAccuracies <- getAccuraciesPerTumorTypeSize(predictionsMM = predictionsMinority,
                                                        metaDataRef = minority$metaDataRef,
                                                        classColumn = higherClassColumn,
                                                        probabilityThreshold = probabilityThresholdMinority)

    majorityAccuracies <- getAccuraciesPerTumorTypeSize(predictionsMM = predictionsMajority,
                                                        metaDataRef = majority$metaDataRef,
                                                        classColumn = higherClassColumn,
                                                        probabilityThreshold = probabilityThresholdMajority)
  } else {
    minorityAccuracies <- getAccuraciesPerTumorTypeSize(predictionsMM = predictionsMinority,
                                                        metaDataRef = minority$metaDataRef,
                                                        classColumn = classColumn,
                                                        probabilityThreshold = probabilityThresholdMinority)

    majorityAccuracies <- getAccuraciesPerTumorTypeSize(predictionsMM = predictionsMajority,
                                                        metaDataRef = majority$metaDataRef,
                                                        classColumn = classColumn,
                                                        probabilityThreshold = probabilityThresholdMajority)

  }

  minorityAccuracies$classifier <- "Minority Classifier"

  majorityAccuracies$classifier <- "Majority Classifier"

  fractionsCorrectTotal <- base::rbind(minorityAccuracies, majorityAccuracies)

  return(fractionsCorrectTotal)
}
