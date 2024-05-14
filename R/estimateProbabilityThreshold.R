#' Calculate metrics to estimate probability score threshold
#'
#' This function calculates the sensitivity, specificity, precision and recall at
#' different probability score thresholds for sample classification.
#' A sample is classified as positive when it receives a classification (classification probability > threshold).
#' On the contrary, a sample is classified as negative when it does not receive a classification
#' (prediction probability < threshold). A true positive is a classified sample for which the classification is
#' correct, a true negative a non-classified sample that had a wrong classification
#'
#' @param predictionsMM  Dataframe showing the top 3 classifications for the tumor (sub)type,
#' together with their probability scores and the original pathology label (originalCall)
#' @param returnPlot Do you want to obtain the data or (FALSE) or get the resulting plot (TRUE)?
#' @param interceptionPoint At which threshold would you like to draw the line that intersects the sensitivity and specificity lines?
#' @return If returnPlot = FALSE: Dataframe containing the number of true positives ($TP), false positives ($FP),
#' true negatives ($TN) and false negatives ($FN) at different probability score thresholds.
#'  Calculated from the metrics named above are the sensitivity ($sensitivity),
#'  specificity ($specificity), precision ($precision) and recall ($recall).
#'
#' If returnPlot = TRUE: Plot showing the sensitivity and specificity at the different probability score thresholds.
#' @export
#'
#'

estimateProbabilityThreshold <- function(predictionsMM,
                                         returnPlot,
                                         interceptionPoint) {

  for (cutoff in base::seq(from = 0, to = 1.01, by = 0.01)) {
    predictionsMMFiltered <- predictionsMM %>% dplyr::filter(probability1 > cutoff)

    TP <- predictionsMMFiltered %>% dplyr::filter(originalCall == predict) %>% base::nrow()
    FP <- predictionsMMFiltered %>% dplyr::filter(originalCall != predict) %>% base::nrow()
    TN <- predictionsMM %>% dplyr::filter(probability1 <= cutoff,
                                   originalCall != predict) %>% base::nrow()
    FN <- predictionsMM %>% dplyr::filter(probability1 <= cutoff,
                                          originalCall == predict) %>% base::nrow()

    #recall <- nrow(predictionsMMFiltered) / nrow(predictionsMM)
    binaryScoreDF <- base::data.frame(cutoff = cutoff,
                                TP = TP,
                                FP = FP,
                                TN = TN,
                                FN = FN)
    if (cutoff == 0) {
      totalbinaryScoreDF <- binaryScoreDF
    } else {
      totalbinaryScoreDF <- base::rbind(totalbinaryScoreDF,
                                  binaryScoreDF)
    }
  }
  totalbinaryScoreDF$sensitivity <- totalbinaryScoreDF$TP/ (totalbinaryScoreDF$TP + totalbinaryScoreDF$FN)

  totalbinaryScoreDF$specificity <- totalbinaryScoreDF$TN/ (totalbinaryScoreDF$TN + totalbinaryScoreDF$FP)

  totalbinaryScoreDF$precision <- totalbinaryScoreDF$TP / (totalbinaryScoreDF$TP + totalbinaryScoreDF$FP)

  totalbinaryScoreDF$recall <- (totalbinaryScoreDF$TP + totalbinaryScoreDF$FP) /
    (totalbinaryScoreDF$TP + totalbinaryScoreDF$FP + totalbinaryScoreDF$TN + totalbinaryScoreDF$FN)

  totalbinaryScoreDF$precision[base::is.na(totalbinaryScoreDF$precision)] <- 1

  if (returnPlot == T) {
    totalbinaryScoreDFGGPlot <- totalbinaryScoreDF %>% tidyr::pivot_longer(cols = sensitivity:recall,
                                                            values_to = "counts")
    thresholdPlot <- totalbinaryScoreDFGGPlot %>% dplyr::filter(name %in% c("sensitivity",
                                                                 "specificity")) %>%
      ggplot2::ggplot(
        ggplot2::aes(x = cutoff,
            linetype = name,
            y = counts
        )) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = interceptionPoint,
                 lty = 1,
                 color = "red") +
      ggplot2::theme_classic() +
      ggplot2::ylab("Sensitivity / Specificity") +
      ggplot2::xlab("Threshold") +
      ggplot2::theme(legend.title= ggplot2::element_blank(),
            plot.margin =  ggplot2::unit(c(0.8,0.8,0.8,0.8), "cm")) +
      ggplot2::scale_linetype_discrete(name = "name", labels = c("Sensitivity","Specificity")) +

      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02)))
    return(thresholdPlot)

  } else {
    return(totalbinaryScoreDF)
  }

}
