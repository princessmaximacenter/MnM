#' Title
#'
#' @param fractionsCorrectTotal
#'
#' @return plot with the separate classifier accuracies per frequency bin.
#' @export
#' @import ungeviz
#'
plotSeparateClassifierAccuracies <- function(fractionsCorrectTotal) {

  if (require("ungeviz") == F) {
    devtools::install_github("wilkelab/ungeviz")
  }
  fractionsCorrectTotal %>%
      ggplot(aes(
      x = nCases,
      y = fractionCorrect,
      col = classifier,
      #fill = classifier
    )) +
    #geom_point(shape = 15) +
    geom_hpline(stat = "identity", size = 0.5) +
    #geom_point() +
    # geom_bar(stat = "identity",
    #         position = "dodge",
    #         color = "black",
    #         width = 0.6
    #         ) +
    theme_classic() +
    ylim(0, 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.05, hjust=1)) +
    labs(x = "Number of patients per tumor type (n)",
         y = "Accuracy") +
    theme(legend.title = element_blank()) +
    theme(#axis.text.x = element_text(hjust = -0.2),
      axis.title.x = element_text(vjust = -1.8),
      axis.title.y = element_text(vjust = 2),
      legend.position = "none",
      plot.margin = unit(c(0.2, 0.6, 0.5, 0.5), "lines")) +
    scale_color_manual(values = c("Minority Classifier" = "#f8766d",
                                  "Majority Classifier" =  "#00bfc4"))
# e78784
  # 7daff9
}
