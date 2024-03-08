#' Calculate performance of M&M on total train and test set
#'
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the reference cohort.
#' @param metaDataTest Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the test set.
#' @param throwOut Are there samples you would like to remove from the test set due to poor data quality? If so, add their rownames here.
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' @param crossValidation Specify whether the results are from the cross-validation setup or not.
#' @param nModels  How many models were created for the majority voting system?
#' @param nSeeds How many seeds was the cross-validation setup run with?
#' @param rounding Do you want rounded numbers for the performance scores?
#' @param probabilityThreshold What is the probability score threshold you would like to use to call a classification 'confident'?
#'
#' @return Dataframe containing the average performance of tumor classifications within a certain frequency range (nCases).
#' Included are the averages for the percentage of correctly and incorrectly classified samples ($meanFractionCorrect and $meanFractionIncorrect),
#' correctly and incorrectly classified 'confident' samples ($meanFractionCorrectFiltered and $meanFractionIncorrectFiltered),
#' averaged precision values for the average of all tumor entities within the frequency range ($meanPrecision),
#' averaged recall values for the average of all tumor entities within the frequency range ($meanRecall),
#' and averaged F1 scores for the average of all tumor entities within the frequency range ($meanF1).
#' Please note that the precision, F1 and recall are calculated for the confident sample classifications only.
#'The total amount of samples within each frequency range ($meanSamples) is also specified.
calculateMeanAndSDAccuracy <- function(classColumn,
                                       higherClassColumn,
                                       minorityDir,
                                       majorityDir,
                                       metaDataRef,
                                       metaDataTest,
                                       throwOut,
         subtype = F,
         crossValidation = T,
         nModels,
         nSeeds,
         rounding,
         probabilityThreshold
         ) {

  for (i in seq(1:nSeeds)) {
    if (crossValidation == T & nSeeds > 1) {
      minorityDoc <- paste0(minorityDir, "seed",i, "/crossValidationMinorityResults.rds")
      majorityDoc <- paste0(majorityDir, "seed",i, "/crossValidationMajorityResults.rds")
    } else {
      minorityDoc <- minorityDir
      majorityDoc <- majorityDir
    }

    minority <- readRDS(minorityDoc)
    majority <- readRDS(majorityDoc)

    predictionsMMFinalList <- integrateMM(minority = minority,
                                      majority = majority,
                                      nModels = nModels,
                                      subtype = subtype,
                                      metaDataRef = metaDataRef,
                                      classColumn = classColumn,
                                      higherClassColumn = higherClassColumn,
                                      crossValidation = crossValidation
    )

    predictionsMMFinal <- predictionsMMFinalList$predictionsMMFinal
    if (crossValidation == F ) {
      predictionsMMFinal %<>% filter(rownames(.) %notin% throwOut)
      if (subtype == F) {
      predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), higherClassColumn]
      } else {
      predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), classColumn]
      }

    }

    if (subtype == T) {
      fractionsCorrect <- getAccuraciesPerTumorTypeSize(predictionsMMFinal,
                                                        metaDataRef = metaDataRef,
                                                        classColumn = classColumn,
                                                        rounding = rounding,
                                                        probabilityThreshold = probabilityThreshold)
    } else {
      fractionsCorrect <- getAccuraciesPerTumorTypeSize(predictionsMMFinal,
                                                        metaDataRef = metaDataRef,
                                                        classColumn = higherClassColumn,
                                                        rounding = rounding,
                                                        probabilityThreshold = probabilityThreshold)

    }
    fractionsCorrect$seed <- i
    if (i == 1) {
      accuracyDF <- fractionsCorrect
    } else {
      accuracyDF <- rbind(accuracyDF, fractionsCorrect)
    }
  }

  meanNumbers <- accuracyDF %>% group_by(nCases) %>%
    summarise(
      meanFractionCorrect = mean(fractionCorrect),
      meanFractionCorrectFiltered = mean(fractionCorrectFiltered),
      meanFractionIncorrect = mean(1 - fractionCorrect),
      meanFractionIncorrectFiltered = mean(1 - fractionCorrectFiltered),
      sdFractionCorrect = sd(fractionCorrect),
      sdFractionCorrectFiltered = sd(fractionCorrectFiltered),
      meanCasesFiltered = round(mean(nSamplesFiltered), digits = 0),
      meanPrecision = mean(Precision),
      meanF1 = mean(F1),
      meanFractionCorrect2 = mean(fractionCorrect2),
      meanFractionCorrect3 = mean(fractionCorrect3),
      sdFractionCorrect2 = sd(fractionCorrect2),
      sdFractionCorrect3 = sd(fractionCorrect3),
      meanRecall = mean(Recall),
      medianF1 = median(F1),
      sdPrecision = sd(Precision),
      sdF1 = sd(F1),
      sdRecall = sd(Recall),
      meanSamples = mean(nSamples)
    )

  return(meanNumbers)
}
