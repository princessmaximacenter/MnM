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
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#'
#' @return Dataframe containing the necessary data points to generate PR-plots for M&M versus other classifiers.
#' @export
precisionRecallPlotData <- function(
                                    minorityDir,
                                    majorityDir,
                                    minorityDirTest,
                                    majorityDirTest,
                                    metaDataTest,
         #predictionsMMTest,
    otherClassifierResultsTrain,
    otherClassifierResultsTest,
         otherClassifierName,
         originalCallColumnOtherClassifier,
    otherClassifierPredictionColumn,
    scoreName,
    subtype
         ) {

  allDirsMinority <- list.dirs(minorityDir, recursive = F)
  allDirsMajority <- list.dirs(majorityDir, recursive = F)
  selectedDirsMinority <- allDirsMinority[grep("seed", allDirsMinority)]
  selectedDirsMajority <- allDirsMajority[grep("seed", allDirsMajority)]

  if (length(selectedDirsMinority) != length(selectedDirsMajority)) {
    stop("The number of models for the minority and majority classifier are not the same.
         Please check your models within the minorityDir and majorityDir")
  } else if (!identical(sub(majorityDir, "", selectedDirsMajority), sub(minorityDir, "", selectedDirsMinority)) ) {
    stop(paste("It seems that classifications from the minority and majority classifier have not been run using the same seed.",
               "\nPlease make sure you run the crossvalidation with the same seed for complementary classifications."))
  }

  for (i in seq(1:length(selectedDirsMajority))) {
    minorityDoc <- paste0(selectedDirsMinority[i],"/crossValidationMinorityResults.rds")
    majorityDoc <- paste0(selectedDirsMajority[i],"/crossValidationMajorityResults.rds")
    minority <- readRDS(minorityDoc)
    majority <- readRDS(majorityDoc)

    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          subtype = subtype
    )

    predictionsMMFinal <- predictionsMMFinalList$predictionsMMFinal

    if (otherClassifierName %in% c("CNS methylation classifier",
                                   "Sarcoma classifier"
    )) {
      predictionsMMFinal[predictionsMMFinal$probability1 < 0.3,"predict"] <- "Geen classificatie"

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
    summarise(Sensitivity = mean(Sensitivity),
              Specificity = mean(Specificity),
              Precision = mean(Precision),
              Recall = mean(Recall),
              minPrecision = min(Precision),
              maxPrecision = max(Precision),

    )
  ourDFTrain$trainOrTest <- "Train"
  #colnames(ourDFTrain) <- gsub("mean", "", colnames(ourDFTrain))
  minorityDoc <- paste0(minorityDirTest, "/minorityClassifierResult.rds")
  majorityDoc <- paste0(majorityDirTest, "/majorityClassifierResult.rds")

  minorityTest <- readRDS(minorityDoc)
  majorityTest <- readRDS(majorityDoc)

  classColumn <- minorityTest$metaDataRun$classColumn
  higherClassColumn <- minorityTest$metaDataRun$higherClassColumn
  predictionsMMTestList <- integrateMM(minority = minorityTest,
               majority = majorityTest,
               subtype = subtype)


  predictionsMMTest <- predictionsMMTestList$predictionsMMFinal

  if (subtype == T) {
    predictionsMMTest$originalCall <- metaDataTest[rownames(predictionsMMTest),classColumn]
  } else {
    predictionsMMTest$originalCall <- metaDataTest[rownames(predictionsMMTest),higherClassColumn]
  }

  if (otherClassifierName %in% c("CNS methylation classifier",
                                 "Sarcoma classifier"
  )) {
    predictionsMMTest[predictionsMMTest$probability1 < 0.3, "predict"] <- "Geen classificatie"
      }

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
