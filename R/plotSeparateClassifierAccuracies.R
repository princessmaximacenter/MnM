#' Plot the separate accuracies of the Minority, Majority and M&M classifier
#'
#' @param separateClassifierAccuracies Dataframe coming from the function 'combineSeparateClassifierAccuracies'.
#'
#' @return Plot showing the accuracy for different population frequencies separately for the Minority, Majority and M&M classifier.
#' @export
#'
plotSeparateClassifierAccuracies <- function(separateClassifierAccuracies) {



  plotSeparateScores <- separateClassifierAccuracies %>%

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

  if ("M&M" %in% unique(separateClassifierAccuracies$classifier)){
    plotSeparateScores <- plotSeparateScores +
      scale_color_manual(values = c("M&M" = "#606ca5",
                                    "Minority Classifier" = "#f8766d",
                                    "Majority Classifier" =  "#00bfc4"))
  } else {
    plotSeparateScores <- plotSeparateScores +
      scale_color_manual(values = c("Minority Classifier" = "#f8766d",
                                    "Majority Classifier" =  "#00bfc4"))

  }

    return(plotSeparateScores)

}
