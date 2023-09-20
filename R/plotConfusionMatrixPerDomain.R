#' Plot tumor subtype confusion matrix
#'
#' Function to plot the confusion matrix plot for a single domain,
#' showing the reference-prediction combinations for the tumor subtypes.
#' Using a color coding, it is shown which tumor subtypes belong to the same overheading tumor type.
#' IMPORTANT: The order of the colors will be applied to the tumor types in the order of the
#' tumors within the abbreviationSubtype dataframe. Make sure you specify the
#' order of the tumor types well within the abbreviationSubtype, not only concerning the domains,
#' but also the subsequent tumor types order and tumor subtypes order.
#'
#' @param domain Which domain do you want to plot? This name should be present in the
#' domainColumn of the metadata.
#' @param confusionPlotDFSubtype Dataframe specifying how often certain reference-prediction
#' combinations for the tumor subtypes are present.
#' @param nonAvailableTiles  A dataframe containing the empty tiles for a confusion matrix plot.
#' @param abbreviationSubtype Dataframe containing the links between the tumor subtype,
#' the abbreviation required in the plot, the tumor type and the domain.
#' @param domainCol Which colors should we use for each designated tumor type?
#' @param colorTiles Which domain color do we want to use?
#'
#' @return ggplot object containing the tumor subtype confusion matrix.
#' @export
#' @import ggplot2 magrittr
#'
plotConfusionMatrixPerDomain <- function(domain,
                                         confusionPlotDFSubtype,
                                         abbreviationSubtype,
                                         nonAvailableTiles,
                                         domainCol,
                                         colorTiles = "#012695") {

  DomainDF <- confusionPlotDFSubtype %>% filter(Domain == domain)
  abbreviationsDomain <- abbreviationSubtype %>% filter(Domain == domain,
                                                         abbreviation %in% as.character(unique(c(DomainDF$Reference, DomainDF$Prediction))) )
  domainSubtypes <- c("Not classified", unique(abbreviationsDomain$abbreviation))

  DomainDF$Prediction <- as.character(DomainDF$Prediction) %>% factor(. , levels = domainSubtypes)

  DomainDF$Reference <- as.character(DomainDF$Reference) %>% factor(. , levels = domainSubtypes)
  #  DomainDF$Reference <- factor(DomainDF$Reference,
  #                               levels = domainSubtypes)

  "Nieuwe factor levels moeten gemaakt worden omdat we maar een sub-deel gebruiken van de labels."

  nonAvailableTilesDomain <- nonAvailableTiles %>% filter(Reference %in% domainSubtypes,
                                                          Prediction %in% domainSubtypes)
  nonAvailableTilesDomain$Prediction <- as.character(nonAvailableTilesDomain$Prediction) %>%
    factor(., levels = domainSubtypes)

  nonAvailableTilesDomain$Reference <- as.character(nonAvailableTilesDomain$Reference) %>%
    factor(., levels = domainSubtypes)


  notClassifiedDF <- data.frame(Reference = "Not classified",
                                Prediction = domainSubtypes, #[domainSubtypes != "Not classified"],
                                Freq = 0,
                                Domain = NA)


  notClassifiedDF$TumorType <- NA
  for (i in seq(1:nrow(notClassifiedDF))) {
    if (notClassifiedDF$Prediction[i] %in% abbreviationsDomain[,"abbreviation"]) {
      notClassifiedDF$TumorType[i] <- abbreviationsDomain[abbreviationsDomain[,"abbreviation"] == notClassifiedDF$Prediction[i], "TumorType"]
    }
  }

  namesDomain <- unique(notClassifiedDF$TumorType)

  notClassifiedDF$TumorType <- factor(notClassifiedDF$TumorType, levels = namesDomain)

  confusionPlot <- ggplot(data = DomainDF, aes(y = Prediction,
                              x = Reference,
                              label = Freq)) +

    geom_tile(color = "black",
              fill = colorTiles) + coord_equal() +

    geom_text(color = "white") +
    guides(fill=F) + # removing legend for `fill`
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +

    geom_tile(data = (nonAvailableTilesDomain), fill = "white", color = "lightgrey") +
    scale_y_discrete(drop=FALSE) + scale_x_discrete(drop = FALSE) +
    geom_tile(data = notClassifiedDF, aes(fill = TumorType),
              color = "black") +
    scale_fill_manual(values = c("grey", domainCol),
                      breaks = namesDomain)

  return(confusionPlot)
}
