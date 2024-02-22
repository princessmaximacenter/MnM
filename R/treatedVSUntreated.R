#' Title
#'
#' @param metaDataRef
#' @param predictionsMM
#' @param probabilityScore
#'
#' @return
#' @export
#'
#' @examples
treatedVSUntreated <- function(metaDataRef,
                               predictionsMM,
                               probabilityScore,
                               testSet
                               ) {
  labelTypes <- c("Yes", "No")

  for (i in seq(1:length(labelTypes))) {
    samples <-  metaDataRef %>% filter(SystemicTreatment == labelTypes[i]) %>% rownames()
    errorsSamples <- predictionsMM[samples,] %>% filter(originalCall != predict) %>% nrow()
    correctSamples <- length(samples) - errorsSamples
    samplesFiltered <- predictionsMM[samples,] %>% filter(probability1 > probabilityScore) %>% rownames()
    errorsFiltered <- predictionsMM[samplesFiltered,] %>% filter(originalCall != predict)%>% nrow()
    correctSamplesFiltered <- length(samplesFiltered) - errorsFiltered
    recall <- round(length(samplesFiltered) / length(samples), digits = 2)
    sampleTypeDF <- data.frame(labelType = labelTypes[i],
                               numberSamples = length(samples),
                               correctSamples = correctSamples,
                               errorsSamples = errorsSamples,
                               accuracy = round(correctSamples / length(samples), digits = 3),
                               numberSamplesFiltered = length(samplesFiltered),
                               correctSamplesFiltered = correctSamplesFiltered,
                               errorsFiltered = errorsFiltered,
                               precision = round(correctSamplesFiltered / length(samplesFiltered),
                                                 digits = 3),
                               recall = recall
    )

    if (i == 1) {
      sampleTypeDFTotal <- sampleTypeDF

    } else {
      sampleTypeDFTotal <- rbind(sampleTypeDFTotal,
                                 sampleTypeDF)
    }

  }

  return(sampleTypeDFTotal)
}
