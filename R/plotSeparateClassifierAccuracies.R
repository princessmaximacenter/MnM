#' Plot the separate accuracies of the Minority classifier, Majority classifier and M&M classifier
#'
#' @param separateClassifierAccuracies Dataframe obtained from the function 'combineSeparateClassifierAccuracies'.
#'
#' @return Plot showing the accuracy for different population frequencies separately for the Minority, Majority and M&M classifier.
#' @export
#'
plotSeparateClassifierAccuracies <- function(separateClassifierAccuracies) {


  if (base::requireNamespace("ungeviz") == F) {
    remotes::install_github("fwallis/ungeviz")
  }

  plotSeparateScores <- separateClassifierAccuracies %>%

    ggplot2::ggplot(
      ggplot2::aes(
      x = whichTop,
      y = accuracy,
      col = classifier,
      group = classifier
    )) +
    ungeviz::geom_hpline(
      stat = "identity",
      size = 1) +

    ggplot2::theme_classic() +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(x = "Number of patients per tumor type (n)",
         y = "Accuracy") +

    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(vjust = -1, size = 20),
      axis.title.y = ggplot2::element_text(vjust = 2, size = 20),

      axis.text.x = ggplot2::element_text(size = 15,  hjust=1, vjust = 0.05, angle = 90),
      axis.text.y = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(size = 11),
      #legend.position = "none",
      plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +
    ggplot2::facet_grid( ~ nCases) +
    ggplot2::geom_vline(xintercept = 2.6, linetype = 2)

  if ("M&M" %in% base::unique(separateClassifierAccuracies$classifier)){
    plotSeparateScores <- plotSeparateScores +
      ggplot2::scale_color_manual(values = c("M&M" = "#606ca5",
                                    "Minority Classifier" = "#f8766d",
                                    "Majority Classifier" =  "#00bfc4"))
  } else {
    plotSeparateScores <- plotSeparateScores +
      ggplot2::scale_color_manual(values = c("Minority Classifier" = "#f8766d",
                                    "Majority Classifier" =  "#00bfc4"))

  }

    return(plotSeparateScores)

}
