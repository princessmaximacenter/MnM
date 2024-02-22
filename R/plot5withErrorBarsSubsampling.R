plot5withErrorBarsSubsampling <- function(classColumn,
                                          higherClassColumn,
                                          minorityDir,
                                          majorityDir,
                                          minorityDirTest,
                                          majorityDirTest,
                                          metaDataRef,
                                          metaDataTest = NA,
                                          otherDataSets,
                                          subsamplePercentage,
                                          #   otherDataSetsTest,
                                          throwOut = NA,
                                          crossValidation = T,
                                          nModels,
                                          nSeeds
) {

  predictionsMMAverageList <- combineSeedPredictions(nSeeds = nSeeds,
                                                     minorityDir = minorityDir,
                                                     majorityDir = majorityDir,
                                                     subtype = F,
                                                     classColumn = classColumn,
                                                     higherClassColumn = higherClassColumn,
                                                     metaDataRef = metaDataRef
  )
  predictionsMMAverage <- predictionsMMAverageList$predictionsMMFinal

  predictionsMMAverageSubtypeList <- combineSeedPredictions(nSeeds = nSeeds,
                                                            minorityDir = minorityDir,
                                                            majorityDir = majorityDir,
                                                            subtype = T,
                                                            classColumn = classColumn,
                                                            higherClassColumn = higherClassColumn,
                                                            metaDataRef = metaDataRef
  )

  predictionsMMAverageSubtype <- predictionsMMAverageSubtypeList$predictionsMMFinal



  minority <- readRDS(minorityDirTest)
  majority <- readRDS(majorityDirTest)

  predictionsMMSubtypeTestList <- integrateMM(minority = minority,
                                              majority = majority,
                                              metaDataRef = metaDataRef,
                                              classColumn = classColumn,
                                              higherClassColumn = higherClassColumn,
                                              nModels = nModels,
                                              subtype = T,
                                              crossValidation = F)

  predictionsMMSubtypeTest <- predictionsMMSubtypeTestList$predictionsMMFinal

  predictionsMMSubtypeTest$originalCall <- metaDataTest[rownames(predictionsMMSubtypeTest),classColumn]

  predictionsMMTestList <- integrateMM(minority = minority,
                                       majority = majority,
                                       metaDataRef = metaDataRef,
                                       classColumn = classColumn,
                                       higherClassColumn = higherClassColumn,
                                       nModels = nModels,
                                       subtype = F,
                                       crossValidation = F)

  predictionsMMTest <- predictionsMMTestList$predictionsMMFinal

  predictionsMMTest$originalCall <- metaDataTest[rownames(predictionsMMTest),higherClassColumn]

  #i <- 1

  for (i in seq(1:length(otherDataSets))) {
    dataSetName <- names(otherDataSets)[i]
    probabilityScoreThreshold <- otherDataSets[[i]]$probabilityScoreThreshold
    subtype <- otherDataSets[[i]]$subtype
    otherDataSet <- otherDataSets[[i]]$data
    subsetName <- otherDataSets[[i]]$subsetName
    trainOrTest <- otherDataSets[[i]]$TrainOrTest
    subsampleNumber <- round(subsamplePercentage * nrow(otherDataSet))
    if (trainOrTest == "Train") {
      crossValidation <- T
      if (subtype == T) {
        predictionsMM <- predictionsMMAverageSubtype
      } else {
        predictionsMM <- predictionsMMAverage
      }
      for (subsampleSelection in seq(1:100)) {
        if (subsampleSelection == 1) {
          set.seed(subsampleSelection)
        }
        selectedSamples <- unique(sample(rownames(otherDataSet), size = subsampleNumber))
        selectedSamplesList <- list(selectedSamples)

        if (subsampleSelection == 1) {
          totalSelectedSamplesList <- selectedSamplesList
        } else {
          totalSelectedSamplesList <- c(totalSelectedSamplesList,
                                        selectedSamplesList)
        }

      }


    } else {
      crossValidation <- F
      totalSelectedSamplesList <- list(rownames(otherDataSet))
      if (subtype == T) {
        predictionsMM <- predictionsMMSubtypeTest
      } else {
        predictionsMM <- predictionsMMTest
      }

    }

    for (subsampleSelection in seq(1:length(totalSelectedSamplesList))) {

      predictionsMnM <- predictionsMM %>% filter(rownames(.) %in% totalSelectedSamplesList[[subsampleSelection]])

      if (dataSetName %in% c("ALLCatchR", "ALLCatchRTest")) {
        errorsMnM <- predictionsMnM %>% filter((predict != originalCall)) %>% nrow()
        predictionsOther <- otherDataSet %>% filter(rownames(.) %in% totalSelectedSamplesList[[subsampleSelection]])

        errorsOther <- predictionsOther %>% filter(originalCall != Prediction) %>% nrow()

        filteredOther <- predictionsOther %>% filter(Confidence == "high-confidence")
        correctFiltered <- filteredOther %>% filter(originalCall == Prediction)
      } else {

        errorsMnM <- predictionsMnM %>% filter((predict != originalCall | probability1 < 0.3)) %>% nrow()

        predictionsOther <- otherDataSet %>% filter(rownames(.) %in% totalSelectedSamplesList[[subsampleSelection]])

        errorsOther <- predictionsOther %>% filter(Disease_sub_class != methClass) %>% nrow()
        filteredOther <- predictionsOther %>% filter(methClassScore > 0.84)
        correctFiltered <- filteredOther %>% filter(Disease_sub_class == methClass)
      }

      precisionOther <- nrow(correctFiltered) / nrow(filteredOther)
      recallOther <- nrow(filteredOther) / nrow(predictionsOther)
      accuracyOther <- (nrow(predictionsOther) - errorsOther) / nrow(predictionsOther)

      filteredMnM <- predictionsMnM %>% filter(probability1 > probabilityScoreThreshold)
      correctFiltered <- filteredMnM %>% filter(originalCall == predict)
      precisionMnM <- nrow(correctFiltered) / nrow(filteredMnM)
      accuracyMnM <- (nrow(predictionsMnM) - errorsMnM) / nrow(predictionsMnM)
      recallMnM <- nrow(filteredMnM) / nrow(predictionsMnM)

      tTestDF <- data.frame(type = c("M&M", "Other"),
                            Subset = subsetName,
                            TrainOrTest = trainOrTest,
                            nSamples = nrow(otherDataSet),
                            errors = c(errorsMnM, errorsOther),
                            accuracy = c(accuracyMnM, accuracyOther),
                            precision = c(precisionMnM, precisionOther),
                            recall = c(recallMnM, recallOther),
                            diffRecall = recallMnM - recallOther,
                            diffPrecision = precisionMnM - precisionOther,
                            diffAccuracy = accuracyMnM - accuracyOther
      )


      if (subsampleSelection == 1) {
        tTestDFTotal <- tTestDF
      } else {
        tTestDFTotal <- rbind(tTestDFTotal,
                              tTestDF)
      }
    }

    #tTestDFTotal$recallDiff <- tTestDFTotal$recallMnM - tTestDFTotal$recallOther
    #tTestDFTotal$accuracyDiff <- tTestDFTotal$accuracyMnM - tTestDFTotal$accuracyOther
    #tTestDFTotal$precisionDiff <- tTestDFTotal$precisionMnM - tTestDFTotal$precisionOther

    if (i == 1) {
      tTestDFCombined <- tTestDFTotal

    } else {
      tTestDFCombined <- rbind(tTestDFCombined,
                               tTestDFTotal)
    }
  }
return(tTestDFCombined)

}
