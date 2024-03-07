#' Plot the results per frequency range of tumor entities
#'
#' Note: This function has been replaced with the more
#' @param meanAndSDPlotTrain Dataframe resulting from function calculateMeanAndSDAccuracy, for training data
#' @param meanAndSDPlotTest Dataframe resulting from function calculateMeanAndSDAccuracy, for test data
#'
#' @return Plot showing the accuracy for train and test set for the different frequencies
#'
plotMeanSDAccuracy <- function(meanAndSDPlotTrain,
                               meanAndSDPlotTest = NA
                               ) {
  meanAndSDPlotTrain$type <- "Train"

  if (!is.na(meanAndSDPlotTest)[1]) {
  meanAndSDPlotTest$type <- "Test"

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
    ylab("Fraction correct samples") +
    geom_bar(
      width = 1,
      stat = "identity",
      col = "black",
      position = "stack",
      fill = "#8bc644") +
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
              alpha = 1) +
    geom_text(#data = totalMeanThing,
      aes(label = paste0("N = ", meanCasesFiltered), y = 0.05), size = 4,
      angle = 90,
      alpha = 1) +
    theme(plot.title = element_text(hjust = 0.5)) +

    labs(x = "Number of patients per tumor type (n)") +
    theme(legend.title = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(vjust = -1.8, size = 25),
          axis.title.y = element_text(vjust = 2, size = 25),
          axis.text = element_text(size = 13, vjust = 0.05),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.spacing = unit(1.2,"lines"),
          plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0))) +
    facet_grid(~nCases, scales = "free_x", space = "free_x")
}
