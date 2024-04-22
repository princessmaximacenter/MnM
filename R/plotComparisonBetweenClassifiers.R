#' Plot comparison between classifiers
#'
#' @param comparisonDF Dataframe containing the accuracy ($meanAccuracy), precision ($meanPrecision),
#' recall ($meanRecall), and their standard deviations for M&M and other classifirs on a subset of the data ($Subset)
#' with a specific number of samples ($numberSamples), either within the training data or test set ($TrainOrTest).
#' @param subsampling Did you use subsampling to generate your comparisonDF?
#'
#' @return plot showing the accuracy, precision and recall for M&M and other classifiers for specific subsets of the dataset will be returned.
#' Included are error bars for the results based on either different seeds used for the cross-valiation results or subsampling.
#'
plotComparisonBetweenClassifiers <- function(comparisonDF, subsampling) {
  if (subsampling == T) {
  comparisonDF <- comparisonDF %>% group_by(Subset, TrainOrTest, type) %>%
    summarise(
      numberSamples = mean(nSamples),
      meanAccuracy = mean(accuracy),
      meanPrecision = mean(precision),
      meanRecall = mean(recall),
      sdAccuracy = sd(accuracy),
      sdPrecision = sd(precision),
      sdRecall = sd(recall)
    )
  }


  comparisonDFLonger <- pivot_longer(comparisonDF, cols = c(meanAccuracy,
                                                                      meanPrecision,
                                                                      meanRecall
  ), names_to = "measurementType",
  values_to = "value")
  #%>% select(all_of(Subset, TrainOrTest, measurementType))

  comparisonDFLonger <- comparisonDFLonger %>% pivot_longer(cols = c("sdAccuracy", "sdPrecision", "sdRecall"),
                                                                   names_to = "whichSD",
                                                                   values_to = "standardDeviation")

  comparisonDFLonger$measurementType <- gsub("mean", "", comparisonDFLonger$measurementType)

  comparisonDFLonger$whichSD <- gsub("sd", "", comparisonDFLonger$whichSD)


  comparisonDFLonger$TrainOrTest <- factor(comparisonDFLonger$TrainOrTest,
                                               levels = c("Train", "Test"))
  comparisonDFLonger$valuePercent <- paste0(round(comparisonDFLonger$value *100, digits = 0), "%")

  comparisonDFLonger <- comparisonDFLonger %>% filter(whichSD == measurementType)


  ggplot(comparisonDFLonger,
         aes(x = TrainOrTest,
             y = value,
             fill = type,
             #color = TrainOrTest,
             label = valuePercent,
             group = type)) +
    theme_classic() +
    geom_errorbar(#data = tTestDFCombinedLonger3,
      aes(ymin= value - standardDeviation, ymax= value + standardDeviation), width=.2,
                  position=position_dodge(.9),
                  col = "black"

    ) +
    geom_col(
      position = position_dodge(0.9),
      color = "black"
    ) +
    scale_fill_manual(values = c("M&M" = "#606ca5",
                                 "Other" = "#bdc6e5")) +
    facet_grid(measurementType  ~ Subset) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       limits = c(0,1.05)) +
    geom_text(position = position_dodge(0.9),
              aes(y = 0.5),
              size = 2.8,
              color = "white"
    ) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          #axis.text = element_text(size = 13),
          axis.text.x = element_text(size = 15,  hjust=1, vjust = 0.05),
          axis.text.y = element_text(size = 15),
          legend.position = "none",
          strip.text = element_text(
            size = 15)

    ) +
    geom_text(data = (comparisonDFLonger %>% filter(type == "M&M",
                                                        measurementType == "Accuracy")),
              aes(label = paste0("N = ", numberSamples), y = 1.04),
              size = 3.5,
              color = "black")

}
