precisionRecallPlotData <- function(#predictionsMMTrain,
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
                                          nModels = nModels,
                                          subtype = subtype,
                                          metaDataRef = metaDataRef,
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
              # minSensitivity = min(Sensitivity),
              #maxSensitivity = max(Sensitivity),
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
