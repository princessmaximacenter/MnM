#' Plot accuracy within the top 3 highest scoring classification labels M&M
#'
#' @param meanAndSDPlotTrain Dataframe containing the average performance of tumor classifications within a certain frequency range (nCases) for the reference cohort.
#' Object obtained from running the function 'calculateMeanAndSDAccuracy'
#' @param meanAndSDPlotTest  Dataframe containing the average performance of tumor classifications within a certain frequency range (nCases) for the test cohort.
#' Object obtained from running the function 'calculateMeanAndSDAccuracy'

#' @return Plot showing the accuracy if not only the top scoring classification label, but also the second and third highest labels are taken into consideration.
#' Plot is segregated into the different class frequencies.
#' @export
#'
plotTop3Accuracy <- function(meanAndSDPlotTrain,
                             meanAndSDPlotTest

                             ) {

  meanAndSDPlotTrain$trainOrTest <- "Reference cohort"
  meanAndSDPlotTest$trainOrTest <- "Test cohort"
  subtype <- meanAndSDPlotTrain$subtype[1]
  top23 <- base::rbind(meanAndSDPlotTrain, meanAndSDPlotTest)
  top23$trainOrTest <- base::factor(top23$trainOrTest, levels = c("Reference cohort", "Test cohort"))

  top23Longer <- top23 %>% tidyr::pivot_longer(cols = c("meanFractionCorrect",
                                                 "meanFractionCorrect2",
                                                 "meanFractionCorrect3"
  ),
  names_to = "topN",
  values_to = "fractionsCorrect"
  )
  top23Longer$fractionCorrectPercent <- base::paste0(base::round(top23Longer$fractionsCorrect * 100, 0), "%")

  top23Longer <- top23Longer %>% tidyr::pivot_longer(cols = c("sdFractionCorrect", "sdFractionCorrect2", "sdFractionCorrect3"),
                                             names_to = "whichFraction",
                                             values_to = "standardDeviation")


  top23Longer$whichFraction <- base::gsub("sd", "mean", top23Longer$whichFraction)
  top23Longer %<>% dplyr::filter(topN == whichFraction)
  top23Longer$topN[base::which(top23Longer$topN == "meanFractionCorrect")] <- "Top 1"
  top23Longer$topN <- base::gsub("meanFractionCorrect", "Top ", top23Longer$topN)


  top3Plot <- ggplot2::ggplot(top23Longer,
                              ggplot2::aes(x = topN,
             y = fractionsCorrect,
             alpha = topN,
             label = fractionCorrectPercent,
             group = topN)) +
    ggplot2::theme_classic() +

    ggplot2::geom_errorbar(
      ggplot2::aes(ymin= fractionsCorrect - standardDeviation, ymax= fractionsCorrect + standardDeviation), width=.2,
      alpha = 1,
                  position=ggplot2::position_dodge(.9),
                  col = "black"

    ) +
    ggplot2::geom_col(position = "dodge",
             col = "black") +
    ggplot2::facet_grid(trainOrTest~nCases, space = "free", scales = "free" ) +
    ggplot2::theme(
          legend.position = "none",
          axis.title.x = ggplot2::element_text(vjust = -1, size = 15),
          axis.title.y = ggplot2::element_text(vjust = 2, size = 15),

          axis.text.x = ggplot2::element_text(size = 10,  hjust=1, vjust = 0.05, angle = 90),
          axis.text.y = ggplot2::element_text(size = 10),
          strip.text = ggplot2::element_text(size = 11),
          plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +
    ggplot2::geom_text(
      ggplot2::aes(y = fractionsCorrect + 0.05),
              position = ggplot2::position_dodge(0.9),
              size = 2.5,
              alpha = 1,
              color = "black") +

    ggplot2::labs(y = "Fraction correct classifications",
         x = "Number of patients per tumor type (n)"
         ) +
    ggplot2::scale_fill_discrete(name = "topN", labels = c("Top 1 classification",
                                                  "Top 2 classifications",
                                                  "Top 3 classifications")) +

    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05)),
                       limits = c(0,1.05))

  if (subtype == T) {
    top3Plot <- top3Plot + ggplot2::labs(x = "Number of patients per tumor subtype (n)")

  }
  return(top3Plot)
}
