#' Calculate performance stratified on recall rate
#'
#' @param predictionsMM Dataframe containing the top 3 final classification labels ($predict{2,3}) with their accompanying probability scores
#' ($probability{1,2,3} and the original diagnosis label ($originalCall).
#' @param probabilityThreshold  What is the probability score threshold you would like to use to call a
#' classification 'confident' for a M&M prediction?
#'
#' @return Dataframe containing the accuracy of tumor entities stratified based on their recall rate.
howManySamplesPassThreshold <- function(predictionsMM,
                                        probabilityThreshold) {

    sequenceHigher <- c(0.76, 0.51, 0.26, 0.01, 0)
    sequenceLower <- c(1.0, 0.75, 0.5, 0.25, 0)

  allTumorNumbers <- predictionsMM %>% select(originalCall) %>% table()

  predictionsMMFiltered <- predictionsMM %>% filter(probability1 > probabilityThreshold)

  predictionsMMFiltered$originalCall <- factor(predictionsMMFiltered$originalCall,
                                                   levels = sort(unique(predictionsMM$originalCall)))

  filteredTumorNumbers <- predictionsMMFiltered %>% select(originalCall) %>% table()

  percentagesThreshold <- filteredTumorNumbers / allTumorNumbers
  for (i in seq(1:length(sequenceHigher))) {

    included <- percentagesThreshold[percentagesThreshold <= sequenceLower[i] & percentagesThreshold >= sequenceHigher[i]]  %>% names()
    correct <- predictionsMM %>% filter(originalCall %in% included,
                                    originalCall == predict) %>% nrow()
    total <- predictionsMM %>% filter(originalCall %in% included) %>% nrow()

    includedDF <- data.frame(samplesIncluded = paste0(sequenceHigher[i], "-", sequenceLower[i], "%"),
                             accuracy = round(correct / total, digits = 2))
    if (i == 1) {
      includedDFTotal <- includedDF

    } else {
      includedDFTotal <- rbind(includedDFTotal,
                               includedDF)
    }
  }
  return(includedDFTotal)
}
