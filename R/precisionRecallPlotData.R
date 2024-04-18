#' Generate data for PR-plot comparison M&M v.s. others
#'
#' @param minorityDir Directory in which the minority model(s) are stored for the cross-validation setup.
#' @param majorityDir Directory in which the majority model(s) are stored for the cross-validation setup.
#' @param predictionsMMTest Classifications of M&M algorithm on the test set, obtained from the function integrateMM.
#' @param otherClassifierResultsTrain Results for other classifier within the reference cohort.
#' @param otherClassifierResultsTest Results for other classifier within the test set.
#' @param otherClassifierName What is the name of the other classifier?
#' @param originalCallColumnOtherClassifier What column in the results for the other classifier contains the original diagnoses for the samples?
#' @param otherClassifierPredictionColumn What column in the results for the other classifier contains the classifications?
#' @param scoreName What column in the results for the other classifier contains the probability scores?
#' @param nSeeds How many seeds was the cross-validation setup run with?
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the reference cohort.
#' @param nModels How many models should be created for the majority voting system?
#' @param classColumn  Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn  Column in the metadata file that contains the tumor type labels.
#'
#' @return Dataframe containing the necessary data points to generate PR-plots for M&M versus other classifiers.
#' @export
precisionRecallPlotData <- function(
                                    minorityDir,
                                    majorityDir,
         predictionsMMTest,
    otherClassifierResultsTrain,
    otherClassifierResultsTest,
         otherClassifierName,
         originalCallColumnOtherClassifier,
    otherClassifierPredictionColumn,
    scoreName,
    nSeeds = 10,
    subtype,
    metaDataRef,
    nModels,
    classColumn,
    higherClassColumn
         ) {

  for (i in seq(1:nSeeds)) {
      minorityDoc <- paste0(minorityDir, "seed",i, "/crossValidationMinorityResults.rds")
      majorityDoc <- paste0(majorityDir, "seed",i, "/crossValidationMajorityResults.rds")

    minority <- readRDS(minorityDoc)
    majority <- readRDS(majorityDoc)

    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          subtype = subtype,
                                          classColumn = classColumn,
                                          higherClassColumn = higherClassColumn,
                                          crossValidation = T
    )

    predictionsMMFinal <- predictionsMMFinalList$predictionsMMFinal

    if (otherClassifierName %in% c("CNS methylation classifier",
                                   "Sarcoma classifier"
    )) {
      predictionsMMFinal[predictionsMMFinal$probability1 < 0.3,"predict"] <- "Geen classificatie"
      predictionsMMTest[predictionsMMTest$probability1 < 0.3, "predict"] <- "Geen classificatie"
    }
    if (i == 1) {

      ourDFTrain <- PRPlotComparison(predictionsMM = predictionsMMFinal,
                                     otherClassifierResults = otherClassifierResultsTrain,
                                     otherClassifierName = otherClassifierName,
                                     originalCallColumnOtherClassifier = originalCallColumnOtherClassifier,
                                     otherClassifierPredictionColumn = otherClassifierPredictionColumn,
                                     scoreName = scoreName
      )

      ourDFTrain$seed <- i
      ourDFTrainTotal <- ourDFTrain

    } else {
      otherClassifierResultsMnM <- predictionsMMFinal %>% filter(rownames(.) %in% rownames(otherClassifierResultsTrain))
      ourDFTrain <- estimateProbabilityThreshold(otherClassifierResultsMnM)
      ourDFTrain$type <- "M&M"
      colnames(ourDFTrain)[6:10] <- c("Sensitivity", "Specificity",
                                      "Precision", "Recall", "type")

      ourDFTrain$seed <- i
      ourDFTrainTotal <- rbind(ourDFTrainTotal,
                               ourDFTrain)
    }
  }

  ourDFTrain <- ourDFTrainTotal %>% group_by(cutoff, type) %>%
    summarise(meanSensitivity = mean(Sensitivity),
              meanSpecificity = mean(Specificity),
              meanPrecision = mean(Precision),
              meanRecall = mean(Recall),
              minPrecision = min(Precision),
              maxPrecision = max(Precision),

    )
  ourDFTrain$trainOrTest <- "Train"
  colnames(ourDFTrain) <- gsub("mean", "", colnames(ourDFTrain))
  ourDFTest <- PRPlotComparison(predictionsMM = predictionsMMTest,
                                otherClassifierResults = otherClassifierResultsTest,
                                otherClassifierName = otherClassifierName,
                                originalCallColumnOtherClassifier = originalCallColumnOtherClassifier,
                                otherClassifierPredictionColumn = otherClassifierPredictionColumn,
                                scoreName = scoreName
  )
  ourDFTest$minPrecision <- ourDFTest$Precision
  ourDFTest$maxPrecision <- ourDFTest$Precision


  ourDFTest$trainOrTest <- "Test"
  dataPR <- rbind(ourDFTrain,
                 ourDFTest[,colnames(ourDFTrain)])
  dataPR$trainOrTest <- factor(dataPR$trainOrTest, levels = c("Train", "Test"))
  dataPR$type <- factor(dataPR$type, levels = c("M&M", otherClassifierName))
  return(dataPR)
}
