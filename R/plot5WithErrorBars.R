#' Plot comparison between classifiers on total results
#'
#' With this function, the overall results will be displayed for the comparison between different classifiers.
#' Error bars will be generated based on different cross-validation runs.
#'
#' @param minorityDir Directory in which the minority model(s) are stored for the cross-validation setup.
#' @param majorityDir Directory in which the majority model(s) are stored for the cross-validation setup.
#' @param minorityDirTest Directory in which the minority model(s) are store for the test setup.
#' @param majorityDirTest Directory in which the majority model(s) are store for the test setup.
#' @param metaDataTest  Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the test set.
#' @param otherDataSets List containing a description of the other datasets, with their actual data ($data),
#'  whether we're looking at the tumor type of subtype ($subtype),
#'  what is the probability score threshold that should be used for the M&M classifications ($probabilityScoreThreshold),
#'  whether we're looking at cross-validation or test set results ($TrainOrTest, should be "Train" or "Test"),
#'  and what is the name that you would like to attach to the subset of the data ($subsetName).
#' @param plotResults Do you want to obtain the data or (FALSE) or get the resulting plot (TRUE)?
#'
#' @return If plotResults == F, a dataframe containing the accuracy ($meanAccuracy), precision ($meanPrecision),
#' recall ($meanRecall), and their standard deviations for M&M and other classifirs on a subset of the data ($Subset)
#' with a specific number of samples ($numberSamples), either within the training data or test set ($TrainOrTest) will be returned.
#' If plotResults == T, a plot showing the accuracy, precision and recall for M&M and other classifiers for specific subsets of the dataset will be returned.
#' Included are error bars for the results of the 10 seeds used for the cross-valiation results.
#'
#'
plot5WithErrorBars <- function(minorityDir,
                               majorityDir,
                               minorityDirTest,
                               majorityDirTest,
                               metaDataTest = NA,
                               otherDataSets,
                            plotResults = F
) {

  for (i in seq(1:length(otherDataSets))) {

    if (i == 1) {
      if (!is.na(metaDataTest)[1]) {
        print("Checking the performance for the test set based on values provided in dataframe 'metaDataTest'.")

        if (classColumn %notin% colnames(metaDataTest)) {
          print("Please note that the wanted column for the tumor subtype labels cannot be found within 'metaDataTest'.")
          print("Either change the column with the tumor subtype labels to the name: ", classColumn)
          stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
        } else if (higherClassColumn %notin% colnames(metaDataTest)) {

          print("Please note that the wanted column for the tumor type labels cannot be found within 'metaDataTest'.")
          print("Either change the column with the tumor type labels to the name: ", higherClassColumn)
          stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
        } else {
          print(paste0("Found columns ", classColumn, " and ", higherClassColumn, " within metaDataTest specifying the tumor subtype, and tumor type."))

          print("No original call found, adding it from metaDataTest")
        }

      }

    }
    dataSetName <- names(otherDataSets)[i]
    probabilityThreshold <- otherDataSets[[i]]$probabilityScoreThreshold
    subtype <- otherDataSets[[i]]$subtype
    otherDataSet <- otherDataSets[[i]]$data
    subsetName <- otherDataSets[[i]]$subsetName
    trainOrTest <- otherDataSets[[i]]$TrainOrTest
    if (trainOrTest == "Train") {

      crossValidation <- T
      minorityDirectory <- minorityDir
      majorityDirectory <- majorityDir
    } else {
      crossValidation <- F
      minorityDirectory <- minorityDirTest
      majorityDirectory <- majorityDirTest
    }
    #subset <- rownames(otherDataSet)

  meanAndSDPlotTrain <- calculateMeanAndSDOverall(#classColumn,
                                                  #higherClassColumn,
                                                  minorityDir = minorityDirectory,
                                                  majorityDir = majorityDirectory,
                                                  metaDataTest = metaDataTest,
                                                  subtype = subtype,
                                                  subset = rownames(otherDataSet),
                                                  probabilityThreshold = probabilityThreshold
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
