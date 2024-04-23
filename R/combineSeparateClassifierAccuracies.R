#' Get results for separate classifiers as compared to M&M
#'
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param probabilityThreshold What is the probability score threshold you would like to use to call a
#' classification 'confident' for a M&M prediction?
#' @param probabilityThresholdMajority What is the probability score threshold you would like to use to call a
#' classification 'confident' for a Majority Classifier prediction? Default value established based on ROC-curve analyses.
#' @param probabilityThresholdMinority What is the probability score threshold you would like to use to call a
#' classification 'confident' for a Minority Classifier prediction? Default value established based on ROC-curve analyses.
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' @param returnPlot Do you want to obtain the data or (FALSE) or get the resulting plot (TRUE)?
#'
#' @export
#' @return If returnPlot = T: Plot with the separate classifier accuracy (Minority Classifier,
#' Majority Classifier and M&M) stratified on population frequency.
#' If returnPlot = F: Dataframe containing the fraction correct labels within the top1, top2 and top3 classification labels.
#'
combineSeparateClassifierAccuracies <- function(minorityDir,
                                             majorityDir,
                                             probabilityThreshold = 0.8,
                                             probabilityThresholdMajority = 0.9,
                                             probabilityThresholdMinority = 0.73,
                                             subtype,
                                             returnPlot = T
                                             ) {

  if (require("ungeviz") == F) {
    remotes::install_github("fwallis/ungeviz")
  }



  allDirsMinority <- list.dirs(minorityDir, recursive = F)
  allDirsMajority <- list.dirs(majorityDir, recursive = F)
  selectedDirsMinority <- allDirsMinority[grep("seed", allDirsMinority)]
  selectedDirsMajority <- allDirsMajority[grep("seed", allDirsMajority)]

  if (length(selectedDirsMinority) != length(selectedDirsMajority)) {
    stop("The number of models for the minority and majority classifier are not the same.
         Please check your models within the minorityDir and majorityDir that the
         same seeds have been used for the generation of a minority and a majority classifier.")
  } else if (!all.equal(selectedDirsMajority, selectedDirsMinority) ) {
    stop("Please make sure you run the crossvalidation with the same seed for complementary classifications,
         and store them in the same directory.")
  }

  for (i in seq(1:length(selectedDirsMajority))) {
    minorityDoc <- paste0(selectedDirsMinority[i],"/crossValidationMinorityResults.rds")
    majorityDoc <- paste0(selectedDirsMajority[i],"/crossValidationMajorityResults.rds")
    minority <- readRDS(minorityDoc)
    majority <- readRDS(majorityDoc)
    classColumn <- minority$metaDataRun$classColumn
    higherClassColumn <- minority$metaDataRun$higherClassColumn

  fractionsCorrectTotal <- getSeparateClassifierAccuracies(minority = minority,
                                                           majority = majority,
                                                           probabilityThresholdMajority = probabilityThresholdMajority,
                                                           probabilityThresholdMinority = probabilityThresholdMinority,
                                                           subtype = subtype)



  predictionsMMFinalList <- integrateMM(minority = minority,
                                        majority = majority,
                                        subtype = subtype)

  predictionsMMFinal <- predictionsMMFinalList$predictionsMMFinal

  if (subtype == T) {
  fractionsCorrectMM <- getAccuraciesPerTumorTypeSize(predictionsMM =  predictionsMMFinal,
                                                           rounding = F,
                                                           metaDataRef = minority$metaDataRef,
                                                           classColumn = classColumn,
                                                           probabilityThreshold = probabilityThreshold)
  } else {
    fractionsCorrectMM <- getAccuraciesPerTumorTypeSize(predictionsMM =  predictionsMMFinal,
                                                             rounding = F,
                                                             metaDataRef = minority$metaDataRef,
                                                             classColumn = higherClassColumn,
                                                             probabilityThreshold = probabilityThreshold)
  }

  fractionsCorrectMM$classifier <- "M&M"
  fractionsCorrectCombi <- rbind(fractionsCorrectMM, fractionsCorrectTotal)

  colnames(fractionsCorrectCombi) <- sub("fractionCorrect3", "Top3", colnames(fractionsCorrectCombi))
  colnames(fractionsCorrectCombi) <- sub("fractionCorrect2", "Top2", colnames(fractionsCorrectCombi))
  colnames(fractionsCorrectCombi) <- sub("fractionCorrect", "Top1", colnames(fractionsCorrectCombi))

  fractionsCorrectCombi$classifier <- factor(fractionsCorrectCombi$classifier, levels = c("M&M", "Majority Classifier", "Minority Classifier"))
  fractionsCorrectCombi$seed <- i

  if (i == 1) {
    fractionsCorrectTotalCombi <- fractionsCorrectCombi

  } else {
    fractionsCorrectTotalCombi <- rbind(fractionsCorrectTotalCombi,
                                        fractionsCorrectCombi
                                        )
  }
  }

  fractionCorrectTotalPivoted <- fractionsCorrectTotalCombi %>% group_by(nCases, classifier) %>%
    summarise(
      Top1 = mean(Top1),
      Top2 = mean(Top2),
      Top3 = mean(Top3)
    ) %>% pivot_longer(cols = c(Top1, Top3),
                       names_to = "whichTop",
                       values_to = "accuracy")



if (returnPlot == T) {

  plotSeparateScores <- plotSeparateClassifierAccuracies(fractionCorrectTotalPivoted)
  if (subtype == T) {
    plotSeparateScores <- plotSeparateScores + xlab("Number of patients per tumor subtype (n)")
  }

 return(plotSeparateScores)
} else {
  return(fractionCorrectTotalPivoted)
  }
}
