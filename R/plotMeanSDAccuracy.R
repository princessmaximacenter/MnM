#' Plot the results per frequency range of tumor entities
#'
#' Note: This function has been replaced with the more
#' @param meanAndSDPlotTrain Dataframe resulting from function calculateMeanAndSDAccuracy, for reference cohort data
#' @param meanAndSDPlotTest Dataframe resulting from function calculateMeanAndSDAccuracy, for test data
#'
#' @return Plot showing the accuracy for train and test set for the different frequencies
#' @export
#' @import magrittr
plotMeanSDAccuracy <- function(meanAndSDPlotTrain,
                               meanAndSDPlotTest = NA
                               ) {
  meanAndSDPlotTrain %<>% filter(!is.na(meanFractionCorrect))
  meanAndSDPlotTrain$type <- "Train"

  if (!is.na(meanAndSDPlotTest)[1]) {
  meanAndSDPlotTest$type <- "Test"
  meanAndSDPlotTest %<>% filter(!is.na(meanFractionCorrect))

  meanAndSDPlot <- rbind(meanAndSDPlotTrain,
                                meanAndSDPlotTest)


  meanAndSDPlot$type <- factor( meanAndSDPlot$type, levels = c("Train", "Test"))
  } else {
    meanAndSDPlot <- meanAndSDPlotTrain
  }
  meanAndSDPlot$fractionCorrectPercent <- paste0(round(meanAndSDPlot$meanFractionCorrectFiltered * 100, 1), "%")

  ggplot(meanAndSDPlot,
         aes(x = type,
             y = meanFractionCorrectFiltered,
             alpha = type
         )) +
    theme_classic() +
    ylab("Precision") +
    geom_bar(
      width = 1,
      stat = "identity",
      col = "black",
      position = "stack",
      fill = #"#8bc644"
        "#8fb559"
      ) +
    geom_errorbar(
      aes(ymin= meanFractionCorrectFiltered - sdFractionCorrect,
          ymax= meanFractionCorrectFiltered + sdFractionCorrect),
      width=.2,
      position=position_dodge(.9),
      col = "#7F7384"

    ) +
    scale_alpha_manual(values = c("Train" = 1,
                                  "Test" = 0.5)) +
    geom_text(aes(label = fractionCorrectPercent,
                  y =  meanFractionCorrectFiltered - 0.05), size = 4,
             # angle = 90,
              alpha = 1) +
    geom_text(
      aes(label = paste0("N = ", meanCasesFiltered), y = 0.03, hjust = 0), size = 4,
      angle = 90,
      alpha = 1) +

    labs(x = "Number of patients per tumor type (n)") +

  theme(#axis.text.x = element_text(hjust = -0.2),
    axis.title.x = element_text(vjust = -1, size = 20),
    axis.title.y = element_text(vjust = 2, size = 20),

    axis.text.x = element_text(size = 15,  hjust=0.8, vjust = 0.5, angle = 90),
    axis.text.y = element_text(size = 15),
    strip.text = element_text(size = 11),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0))) +
    facet_grid(~nCases, scales = "free_x", space = "free_x")
}
