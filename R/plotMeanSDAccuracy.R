#' Plot the results per frequency range of tumor entities
#'
#' @param meanAndSDPlotTrain Dataframe resulting from function 'calculateMeanAndSDAccuracy', for reference cohort data.
#' @param meanAndSDPlotTest Dataframe resulting from function 'calculateMeanAndSDAccuracy', for test cohort data.
#'
#' @return Plot showing the accuracy for train and test set for the different frequencies.
#' @export
plotMeanSDAccuracy <- function(meanAndSDPlotTrain,
                               meanAndSDPlotTest = NA
                               ) {
  meanAndSDPlotTrain %<>% dplyr::filter(!is.na(meanFractionCorrect))
 subtype <- meanAndSDPlotTrain$subtype[1]

  if (!base::is.na(meanAndSDPlotTest)[1]) {
  meanAndSDPlotTest %<>% dplyr::filter(!is.na(meanFractionCorrect))

  meanAndSDPlot <- base::rbind(meanAndSDPlotTrain,
                                meanAndSDPlotTest)


  meanAndSDPlot$type <- base::factor( meanAndSDPlot$type, levels = c("Train", "Test"))
  } else {
    meanAndSDPlot <- meanAndSDPlotTrain
  }
  meanAndSDPlot$fractionCorrectPercent <- base::paste0(base::round(meanAndSDPlot$meanFractionCorrectFiltered * 100, 1), "%")

  precisionPlot <- ggplot2::ggplot(meanAndSDPlot,
                                   ggplot2::aes(x = type,
             y = meanFractionCorrectFiltered,
             alpha = type
         )) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Precision") +
    ggplot2::geom_bar(
      width = 1,
      stat = "identity",
      col = "black",
      position = "stack",
      fill = #"#8bc644"
        "#8fb559"
      ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin= meanFractionCorrectFiltered - sdFractionCorrect,
          ymax= meanFractionCorrectFiltered + sdFractionCorrect),
      width=.2,
      position=ggplot2::position_dodge(.9),
      col = "#7F7384"

    ) +
    ggplot2::scale_alpha_manual(values = c("Train" = 1,
                                  "Test" = 0.5)) +
    ggplot2::geom_text(
      ggplot2::aes(label = fractionCorrectPercent,
                  y =  meanFractionCorrectFiltered - 0.05), size = 4,
             # angle = 90,
              alpha = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = base::paste0("N = ", meanCasesFiltered), y = 0.03, hjust = 0), size = 4,
      angle = 90,
      alpha = 1) +

    ggplot2::labs(x = "Number of patients per tumor type (n)") +

    ggplot2::theme(#axis.text.x = element_text(hjust = -0.2),
    axis.title.x = ggplot2::element_text(vjust = -1, size = 20),
    axis.title.y = ggplot2::element_text(vjust = 2, size = 20),

    axis.text.x = ggplot2::element_text(size = 15,  hjust=0.8, vjust = 0.5, angle = 90),
    axis.text.y = ggplot2::element_text(size = 15),
    strip.text = ggplot2::element_text(size = 11),
    legend.position = "none",
    plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.1, 0))) +
    ggplot2::facet_grid(~nCases, scales = "free_x", space = "free_x")

  if (subtype == T) {
    precisionPlot <- precisionPlot + ggplot2::labs(x = "Number of patients per tumor subtype (n)")
  }
  return(precisionPlot)
}
