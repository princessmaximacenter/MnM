plot5WithErrorBars <- function(classColumn,
                               higherClassColumn,
                               minorityDir,
                               majorityDir,
                               minorityDirTest,
                               majorityDirTest,
                               metaDataRef,
                               metaDataTest = NA,
                               otherDataSets,
                            #   otherDataSetsTest,
                               throwOut = NA,
                               crossValidation = T,
                               nModels,
                               nSeeds
) {

  for (i in seq(1:length(otherDataSets))) {
    dataSetName <- names(otherDataSets)[i]
    probabilityScoreThreshold <- otherDataSets[[i]]$probabilityScoreThreshold
    subtype <- otherDataSets[[i]]$subtype
    otherDataSet <- otherDataSets[[i]]$data
    subsetName <- otherDataSets[[i]]$subsetName
    trainOrTest <- otherDataSets[[i]]$TrainOrTest
    if (trainOrTest == "Train") {
      crossValidation <- T
    } else {
      crossValidation <- F
    }
    #subset <- rownames(otherDataSet)

  meanAndSDPlotTrain <- calculateMeanAndSDOverall(classColumn = classColumn,
                                                  higherClassColumn = higherClassColumn,
                                                  minorityDir = minorityDir,
                                                  majorityDir = majorityDir,
                                                  minorityDirTest = minorityDirTest,
                                                  majorityDirTest = majorityDirTest,
                                                  metaDataRef = metaDataRef,
                                                  metaDataTest = metaDataTest,
                                                  throwOut = throwOut,
                                                  subtype = subtype,
                                                  crossValidation = crossValidation,
                                                  nModels = nModels,
                                                  nSeeds = nSeeds,
                                                  subset = rownames(otherDataSet),
                                                  probabilityScoreThreshold = probabilityScoreThreshold

                                                  )


  meanAndSDPlotTrain$type <- "M&M"
  meanAndSDPlotTrain$Subset <- subsetName
  meanAndSDPlotTrain$TrainOrTest <- trainOrTest

  totalSamples <- nrow(otherDataSet)
  if (dataSetName %in% c("ALLCatchR", "ALLCatchRTest")) {
    correct <- otherDataSet %>% filter(originalCall == Prediction) %>% nrow()
    accuracy <-  round(correct / totalSamples, digits = 2)
    recalled <- otherDataSet %>% filter(Confidence == "high-confidence") %>% nrow()
    recall <- recalled / totalSamples
    precisioned <- otherDataSet %>% filter(Confidence == "high-confidence",
                                                            originalCall == Prediction) %>% nrow()
    precision <- precisioned / recalled


  } else {

    correct <- otherDataSet %>% filter(!!sym(higherClassColumn) == methClass) %>% nrow()
    accuracy <- round(correct / totalSamples, digits = 2)

    recalled <- otherDataSet %>% filter(methClassScore > 0.84) %>% nrow()
    recall <- recalled / totalSamples
    precisioned <- otherDataSet %>% filter(methClassScore > 0.84,
                                                            !!sym(higherClassColumn) == methClass) %>% nrow()
    precision <- precisioned / recalled

  }

  meanAndSDPlotTrainOtherData <- data.frame(meanAccuracy = accuracy,
                                      meanPrecision = precision,
                                      meanRecall = recall,
                                      sdAccuracy = 0,
                                      sdPrecision = 0,
                                      sdRecall = 0,
                                      type = "Other",
                                      Subset = subsetName,
                                      TrainOrTest = trainOrTest)
  meanAndSDPlotTrainCombi <- rbind(meanAndSDPlotTrain,
                                   meanAndSDPlotTrainOtherData)
  meanAndSDPlotTrainCombi$numberSamples <- totalSamples

 if (i == 1) {
   meanAndSDPlotTotal <- meanAndSDPlotTrainCombi

 } else {
   meanAndSDPlotTotal <- rbind(meanAndSDPlotTotal, meanAndSDPlotTrainCombi)

 }
  }
return(meanAndSDPlotTotal)
}
