#' Plot comparison between classifiers on subsampled results
#'
#' With this function, the results of will be displayed for the comparison between different classifiers.
#' Error bars will be generated based on results for 100 different data subsamples.
#'
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param minorityDir Directory in which the minority model(s) are stored for the cross-validation setup.
#' @param majorityDir Directory in which the majority model(s) are stored for the cross-validation setup.
#' @param minorityDirTest Directory in which the minority model(s) are store for the test setup.
#' @param majorityDirTest Directory in which the majority model(s) are store for the test setup.
#' @param metaDataRef  Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the reference cohort.
#' @param metaDataTest  Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the test set.
#' @param otherDataSets List containing a description of the other datasets, with their actual data ($data),
#'  whether we're looking at the tumor type of subtype ($subtype),
#'  what is the probability score threshold that should be used for the M&M classifications ($probabilityScoreThreshold),
#'  whether we're looking at cross-validation or test set results ($TrainOrTest, should be "Train" or "Test"),
#'  and what is the name that you would like to attach to the subset of the data ($subsetName).
#' @param subsamplePercentage Which percentage of samples would you want to subsample each time?
#' @param nModels  How many models were created for the majority voting system?
#' @param returnPlot Do you want to obtain the data or (FALSE) or get the resulting plot (TRUE)?
#' @param numberSubsamples How often do you want to subsample?
#' @param probabilityScoreMinimum Do you want to say that classifications with a very low probability score should be called 'wrong classifications'?
#' If yes, please specify the probability score this should happen at. This will be used within the results for the methylation classifiers.
#' @return  If returnPlot == F, a dataframe containing the accuracy ($meanAccuracy), precision ($meanPrecision),
#' recall ($meanRecall), and their standard deviations for M&M and other classifirs on a subset of the data ($Subset)
#' with a specific number of samples ($numberSamples), either within the training data or test set ($TrainOrTest) will be returned.
#' If returnPlot == T, a plot showing the accuracy, precision and recall for M&M and other classifiers for specific subsets of the dataset will be returned.
#' Included are error bars for the results of the 10 seeds used for the cross-valiation results.
#'
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
                                          numberSubsamples = 100,
                                          #   otherDataSetsTest,
                                          nModels,
                                          returnPlot = F,
                                          probabilityScoreMinimum = 0.0
) {

  predictionsMMAverageList <- combineSeedPredictions(
                                                     minorityDir = minorityDir,
                                                     majorityDir = majorityDir,
                                                     subtype = F,
                                                     classColumn = classColumn,
                                                     higherClassColumn = higherClassColumn)

  predictionsMMAverage <- predictionsMMAverageList$predictionsMMFinal

  predictionsMMAverageSubtypeList <- combineSeedPredictions(
                                                            minorityDir = minorityDir,
                                                            majorityDir = majorityDir,
                                                            subtype = T,
                                                            classColumn = classColumn,
                                                            higherClassColumn = higherClassColumn)

  predictionsMMAverageSubtype <- predictionsMMAverageSubtypeList$predictionsMMFinal



  minority <- readRDS(paste0(minorityDirTest, "/minorityClassifierResult.rds"))
  majority <- readRDS(paste0(majorityDirTest, "/majorityClassifierResult.rds"))

  predictionsMMSubtypeTestList <- integrateMM(minority = minority,
                                              majority = majority,
                                              classColumn = classColumn,
                                              higherClassColumn = higherClassColumn,
                                              subtype = T)

  predictionsMMSubtypeTest <- predictionsMMSubtypeTestList$predictionsMMFinal

  predictionsMMSubtypeTest$originalCall <- metaDataTest[rownames(predictionsMMSubtypeTest),classColumn]

  predictionsMMTestList <- integrateMM(minority = minority,
                                       majority = majority,
                                       classColumn = classColumn,
                                       higherClassColumn = higherClassColumn,
                                       subtype = F)

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
    #if (trainOrTest == "Train") {
     # crossValidation <- T
      if (subtype == T & trainOrTest == "Train") {
        predictionsMM <- predictionsMMAverageSubtype
      } else if (subtype == T & trainOrTest == "Test") {
        predictionsMM <- predictionsMMSubtypeTest
      } else if (subtype == F & trainOrTest == "Train") {
        predictionsMM <- predictionsMMAverage
      } else if (subtype == F & trainOrTest == "Test") {
        predictionsMM <- predictionsMMTest
      }
      #}
      for (subsampleSelection in seq(1:numberSubsamples)) {
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


    #} else {
     # crossValidation <- F
      #totalSelectedSamplesList <- list(rownames(otherDataSet))
      #if (subtype == T) {
      #  predictionsMM <- predictionsMMSubtypeTest
      #} else {
      #  predictionsMM <- predictionsMMTest
      #}

    #}

    for (subsampleSelection in seq(1:length(totalSelectedSamplesList))) {

      predictionsMnM <- predictionsMM %>% filter(rownames(.) %in% totalSelectedSamplesList[[subsampleSelection]])

      if (dataSetName %in% c("ALLCatchR", "ALLCatchRTest")) {
        errorsMnM <- predictionsMnM %>% filter((predict != originalCall)) %>% nrow()
        predictionsOther <- otherDataSet %>% filter(rownames(.) %in% totalSelectedSamplesList[[subsampleSelection]])

        errorsOther <- predictionsOther %>% filter(originalCall != Prediction) %>% nrow()

        filteredOther <- predictionsOther %>% filter(Confidence == "high-confidence")
        correctFiltered <- filteredOther %>% filter(originalCall == Prediction)
      } else {

        errorsMnM <- predictionsMnM %>% filter((predict != originalCall | probability1 < probabilityScoreMinimum)) %>% nrow()

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

    if (i == 1) {
      tTestDFCombined <- tTestDFTotal

    } else {
      tTestDFCombined <- rbind(tTestDFCombined,
                               tTestDFTotal)
    }
  }
  if (returnPlot == T) {
    return(plotComparisonBetweenClassifiers(tTestDFCombined, subsampling = T))

  } else {
    return(tTestDFCombined)
  }

}
