#' Get results for separate classifiers as compared to M&M
#'
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param probabilityThreshold What is the probability score threshold you would like to use to call a
#' classification 'confident' for a M&M prediction?
#' @param probabilityThresholdMajority What is the probability score threshold you would like to use to call a
#' classification 'confident' for a Majority Classifier prediction?
#' @param probabilityThresholdMinority What is the probability score threshold you would like to use to call a
#' classification 'confident' for a Minority Classifier prediction?
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' @param nModels How many models were created for the majority voting system?
#' @param nSeeds How many seeds was the cross-validation setup run with?
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the reference cohort.
#' @param returnPlot Do you want to obtain the data or (FALSE) or get the resulting plot (TRUE)?
#'
#' @return If returnPlot = T: Plot with the separate classifier accuracy (Minority Classifier,
#' Majority Classifier and M&M) stratified on population frequency.
#' If returnPlot = F: Dataframe containing the fraction correct labels within the top1, top2 and top3 classification labels.
#'
plotSeparateClassifierAccuracies <- function(minorityDir,
                                             majorityDir,
                                             classColumn,
                                             higherClassColumn,
                                             probabilityThreshold,
                                             probabilityThresholdMajority,
                                             probabilityThresholdMinority,
                                             subtype,
                                             nModels,
                                             nSeeds,
                                             metaDataRef,
                                             returnPlot = T
                                             ) {

  if (require("ungeviz") == F) {
    remotes::install_github("fwallis/ungeviz")
  }

  for (i in seq(1:nSeeds)) {
      minorityDoc <- paste0(minorityDir, "seed",i, "/crossValidationMinorityResults.rds")
      majorityDoc <- paste0(majorityDir, "seed",i, "/crossValidationMajorityResults.rds")

    minority <- readRDS(minorityDoc)
    majority <- readRDS(majorityDoc)


  fractionsCorrectTotal <- getSeparateClassifierAccuracies(minority = minority,
                                                           majority = majority,
                                                           crossValidation = T,
                                                           classColumn = classColumn,
                                                           higherClassColumn = higherClassColumn,
                                                           probabilityThresholdMajority = probabilityThresholdMajority,
                                                           probabilityThresholdMinority = probabilityThresholdMinority,
                                                           subtype = subtype,
                                                           nModels = nModels)



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

  if (subtype == T) {
  fractionsCorrectMM <- getAccuraciesPerTumorTypeSize(predictionsMM =  predictionsMMFinal,
                                                           rounding = F,
                                                           metaDataRef = metaDataRef,
                                                           classColumn = classColumn,
                                                           probabilityThreshold = probabilityThreshold)
  } else {
    fractionsCorrectMM <- getAccuraciesPerTumorTypeSize(predictionsMM =  predictionsMMFinal,
                                                             rounding = F,
                                                             metaDataRef = metaDataRef,
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


  plotSeparateScores <- fractionCorrectTotalPivoted %>%

    ggplot(aes(
      x = whichTop,
      y = accuracy,
      col = classifier,
      group = classifier
    )) +
    geom_hpline(
      stat = "identity",
      size = 1) +

    theme_classic() +
    ylim(0, 1) +
    labs(x = "Number of patients per tumor type (n)",
         y = "Accuracy") +

    scale_color_manual(values = c("M&M" = "#606ca5",
                                  "Minority Classifier" = "#f8766d",
                                  "Majority Classifier" =  "#00bfc4"
    )) +

    theme(legend.title = element_blank()) +
    theme(
      axis.title.x = element_text(vjust = -1, size = 20),
      axis.title.y = element_text(vjust = 2, size = 20),

      axis.text.x = element_text(size = 15,  hjust=1, vjust = 0.05, angle = 90),
      axis.text.y = element_text(size = 15),
      strip.text = element_text(size = 11),
      legend.position = "none",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +
    facet_grid( ~ nCases) +
    geom_vline(xintercept = 2.6, linetype = 2)
if (subtype == T) {
  plotSeparateScores <- plotSeparateScores + xlab("Number of patients per tumor subtype (n)")
}
if (returnPlot == T) {
 return(plotSeparateScores)
} else {
  return(fractionCorrectTotalPivoted)
  }
}
