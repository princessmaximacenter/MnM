howManySamplesPassThreshold <- function(predictionsMM,
                                        threshold = 0.8) {

    #sequenceHigher <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0)
    #sequenceLower <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)

    sequenceHigher <- c(0.76, 0.51, 0.26, 0.01, 0)
    sequenceLower <- c(1.0, 0.75, 0.5, 0.25, 0)

  allTumorNumbers <- predictionsMM %>% select(originalCall) %>% table()

  predictionsMMFiltered <- predictionsMM %>% filter(probability1 > threshold)

  predictionsMMFiltered$originalCall <- factor(predictionsMMFiltered$originalCall,
                                                   levels = sort(unique(predictionsMM$originalCall)))

  filteredTumorNumbers <- predictionsMMFiltered %>% select(originalCall) %>% table()

  percentagesThreshold <- filteredTumorNumbers / allTumorNumbers
 # percentages <- c(100, 75, 50, 25, 0)

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
