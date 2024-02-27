#' Plot comparison between classifiers on total results
#'
#' With this function, the overall results will be displayed for the comparison between different classifiers.
#' Error bars will be generated based on different cross-validation runs.
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
#' @param throwOut Are there samples you would like to remove from the test set due to poor data quality? If so, add their rownames here.
#' @param nModels  How many models were created for the majority voting system?
#' @param nSeedsTrain How many seeds was the cross-validation setup run with?
#' @param plotResults Do you want to obtain the data or (FALSE) or get the resulting plot (TRUE)?
#'
#' @return If plotResults == F, a dataframe containing the accuracy ($meanAccuracy), precision ($meanPrecision),
#' recall ($meanRecall), and their standard deviations for M&M and other classifirs on a subset of the data ($Subset)
#' with a specific number of samples ($numberSamples), either within the training data or test set ($TrainOrTest) will be returned.
#' If plotResults == T, a plot showing the accuracy, precision and recall for M&M and other classifiers for specific subsets of the dataset will be returned.
#' Included are error bars for the results of the 10 seeds used for the cross-valiation results.
#'
#'
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
                               nModels,
                            nSeedsTrain,
                            plotResults = F
) {

  for (i in seq(1:length(otherDataSets))) {
    dataSetName <- names(otherDataSets)[i]
    probabilityScoreThreshold <- otherDataSets[[i]]$probabilityScoreThreshold
    subtype <- otherDataSets[[i]]$subtype
    otherDataSet <- otherDataSets[[i]]$data
    subsetName <- otherDataSets[[i]]$subsetName
    trainOrTest <- otherDataSets[[i]]$TrainOrTest
    if (trainOrTest == "Train") {
      nSeeds <- nSeedsTrain
      crossValidation <- T
      minorityDirectory <- minorityDir
      majorityDirectory <- majorityDir
    } else {
      #nSeeds <- 1
      crossValidation <- F
      minorityDirectory <- minorityDirTest
      majorityDirectory <- majorityDirTest
    }
    #subset <- rownames(otherDataSet)

  meanAndSDPlotTrain <- calculateMeanAndSDOverall(classColumn,
                                                  higherClassColumn,
                                                  minorityDir = minorityDirectory,
                                                  majorityDir = majorityDirectory,
                                                  metaDataRef,
                                                  metaDataTest,
                                                  throwOut = throwOut,
                                                  subtype = subtype,
                                                  crossValidation = crossValidation,
                                                  nModels = nModels,
                                                  nSeeds  = nSeeds,
                                                  subset = rownames(otherDataSet),
                                                  probabilityThreshold = probabilityScoreThreshold
                                                  )

  #meanAndSDPlotTrain %<>% filter(nCases == "All")
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
  if (plotResults == T) {
    return(plotComparisonBetweenClassifiers(meanAndSDPlotTotal))

  } else {
    return(meanAndSDPlotTotal)
  }
}
