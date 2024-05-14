#' Plot tumor type or subtype confusion matrix
#'
#' Function to plot the confusion matrix (CM) plot within all domains (domain = NA)
#' or within a single domain (e.g. domain = "Hemato").
#' The CM plots can both be made on the tumor type and subtype level,
#' which depends on the information stored within the 'confusionPlotInfo'.
#' The CM plot shows the reference-classification combinations for the tumor (sub)types.
#' Within tumor subtype CM plots, it is shown which tumor subtypes belong to the same parent tumor type using color coding.

#'
#' @param domain Which domain do you want to plot? This name should be available in the
#' domain column of the abbreviations. If not specified, plot will be made for all tumor domains.
#' @param confusionPlotInfo List resulting from running the function 'getConfusionMatrixPlot'.
#' @param domainColor If one domain is plotted, which colors should we use for each designated tumor type?
#' If not specified, default ggplot colors will be used.
#' @param colorTiles Which domain color tiles do we want to use? Only applicable when argument _domain_ is specified.
#'
#' @return ggplot object containing the tumor type or subtype CM plot.
#' @export
#'
plotConfusionMatrix <- function(domain = NA,
                                         confusionPlotInfo,
                                         domainColor = NA,
                                         colorTiles = "#012695") {

  confusionPlotDF <- confusionPlotInfo$confusionPlotDF
  abbreviations <- confusionPlotInfo$abbreviations
  nonAvailableTiles <- confusionPlotInfo$nonAvailableTiles

  if (base::is.na(domain)) {
    DomainDF <- confusionPlotDF
    nonAvailableTilesDomain <- nonAvailableTiles
  } else {

    DomainDF <- confusionPlotDF %>% dplyr::filter(Domain == domain)
    abbreviationsDomain <- abbreviations %>% dplyr::filter(Domain == domain,
                                                           abbreviation %in% base::as.character(base::unique(c(DomainDF$Reference, DomainDF$Prediction))) )
    abbreviationsExtra <- abbreviations %>% dplyr::filter(Domain != domain,
                                                          abbreviation %in% base::as.character(base::unique(c(DomainDF$Reference, DomainDF$Prediction))),
                                                          Domain != "Not classified")
    domainSubtypes <- c("Not classified", base::unique(abbreviationsDomain$abbreviation), base::unique(abbreviationsExtra$abbreviation))
    #domainSubtypes <- c(unique(abbreviationsDomain$abbreviation))
    DomainDF$Prediction <- base::as.character(DomainDF$Prediction) %>% base::factor(. , levels = domainSubtypes)

    DomainDF$Reference <- base::as.character(DomainDF$Reference) %>% base::factor(. , levels = domainSubtypes)
    #  DomainDF$Reference <- factor(DomainDF$Reference,
    #                               levels = domainSubtypes)

    "Nieuwe factor levels moeten gemaakt worden omdat we maar een sub-deel gebruiken van de labels."

    nonAvailableTilesDomain <- nonAvailableTiles %>% dplyr::filter(Reference %in% domainSubtypes,
                                                                   Prediction %in% domainSubtypes)
    nonAvailableTilesDomain$Prediction <- base::as.character(nonAvailableTilesDomain$Prediction) %>%
      base::factor(., levels = domainSubtypes)

    nonAvailableTilesDomain$Reference <- base::as.character(nonAvailableTilesDomain$Reference) %>%
      base::factor(., levels = domainSubtypes)


    notClassifiedDF <- base::data.frame(
      Prediction = domainSubtypes,
      Reference = "Not classified",
      Freq = 0,
      Domain = NA)


    notClassifiedDF[,"TumorType"] <- NA

    for (i in base::seq(1:base::nrow(notClassifiedDF))) {
      if (notClassifiedDF$Prediction[i] %in% abbreviationsDomain[,"abbreviation"] & confusionPlotInfo$subtype == T) {
        notClassifiedDF$TumorType[i] <- abbreviationsDomain[abbreviationsDomain[,"abbreviation"] == notClassifiedDF$Prediction[i], confusionPlotInfo$higherClassColumn]
      }
    }

    namesDomain <- base::unique(notClassifiedDF$TumorType)

    notClassifiedDF$TumorType <- base::factor(notClassifiedDF$TumorType, levels = namesDomain)
  }

  if (!base::is.na(domain)) {
    confusionPlot <- ggplot2::ggplot(data = DomainDF,
                                     ggplot2::aes(y = Prediction,
                                                          x = Reference,
                                                          label = Freq)) +

      ggplot2::geom_tile(color = "black",
                         fill = colorTiles) + ggplot2::coord_equal()
  } else {
    confusionPlot <- ggplot2::ggplot(data = DomainDF,
                                     ggplot2::aes(y = Prediction,
                                                          x = Reference,
                                                          label = Freq,
                                                          fill= Domain)) +

      ggplot2::geom_tile(color = "black") + ggplot2::coord_equal()

  }
  confusionPlot <- confusionPlot +


    ggplot2::geom_text(color = "white") +
    ggplot2::guides(fill=F) + # removing legend for `fill`
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
                   axis.text = ggplot2::element_text(size = 14),
                   axis.title = ggplot2::element_text(size = 18),
                   plot.margin = ggplot2::unit(c(0.8,0.8,0.8,0.8), "cm")
    ) +

    ggplot2::geom_tile(data = (nonAvailableTilesDomain), fill = "white", color = "lightgrey") +
    ggplot2::scale_y_discrete(drop=FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +

    ggplot2::labs(y = "Classification")

  if (!is.na(domain)) {
    confusionPlot <- confusionPlot +
      ggplot2::geom_tile(data = notClassifiedDF,
                         ggplot2::aes(fill = TumorType),
                         color = "black")
  }

  if (!is.na(domainColor)[1]) {
    confusionPlot <- confusionPlot +
      ggplot2::scale_fill_manual(values = c("grey", domainColor),
                                 breaks = namesDomain)
  }

  return(confusionPlot)
}
