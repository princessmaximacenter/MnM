#' Title
#'
#' @param domain
#' @param ggplotDF
#' @param nonAvailableTiles
#' @param namesDomain
#' @param domainSubtypes
#' @param domainCol
#' @param colorTiles
#'
#' @return
#' @export
#'
plotConfusionMatrixPerDomain <- function(domain,
                                         ggplotDFSubtype,
                                         abbreviationsSubtype,
                                         nonAvailableTiles,
                                         domainCol,
                                         colorTiles = "#012695") {

  abbreviationsDomain <- abbreviationsSubtype %>% filter(Domain == domain,
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
    scale_fill_manual(values = domainCol,
                      breaks = namesDomain)

  return(confusionPlot)
}
