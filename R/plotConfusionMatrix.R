#' Create confusion matrix plot with the three domains.
#'
#' This function can plot a confusion matrix both for the tumor type
#' and the tumor subtype level.
#'
#' @param confusionPlotDF Dataframe specifying how often certain reference-prediction
#' combinations are present.
#' @param nonAvailableTiles A dataframe containing the empty tiles for a confusion matrix plot.
#'
#' @return ggplot object with the confusion matrix of the three domains.
#' @export
#' @import ggplot2
#'
plotConfusionMatrix <- function(confusionPlotDF,
                                nonAvailableTiles) {

  confusionPlot <- ggplot(confusionPlotDF, aes(y = Prediction,
                                               x = Reference, fill= Domain, label = Freq)) +
    geom_tile(color = "black") + coord_equal() +
    geom_text(color = "white") +
    #   scale_fill_continuous(low="white", high="steelblue",
    #                         guide="colorbar",na.value="white") +
    guides(fill=F) + # removing legend for `fill`
    theme_bw() +
    #labs(title = paste(whichClassification, "Classified Samples")) + #Add Title
    # geom_text(aes(label=Freq), color="black") + # Add values to tiles
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = c("Hemato" = "#880808",
                                 "Solid" =  "#D1944A",
                                 "Neuro" = "#012695"),
                      labels = c("Blood tumors", 'Neurological tumors', "Solid tumors")) +
    scale_y_discrete(drop=FALSE) + scale_x_discrete(drop = FALSE) +
    geom_tile(data = nonAvailableTiles, fill = "white", color = "lightgrey")

  return(confusionPlot)
}
