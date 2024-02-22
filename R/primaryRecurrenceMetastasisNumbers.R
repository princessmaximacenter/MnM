primaryRecurrenceMetastasisNumbers <- function(metaDataRef,
         predictionsMM,
         probabilityScore,
         columnOfInterest,
         labelTypes,
         includeAllRow = F
         ) {

  # recurrences <-  metaDataRef %>% filter(Status == "Recurrence") %>% rownames()
  # errorsRecurrence <- predictionsMMAverage[recurrences,] %>% filter(originalCall != predict)%>% nrow()
  # correctRecurrence <- length(recurrences) - errorsRecurrence
  #
  # metastases <- metaDataRef %>% filter(Status == "Metastasis") %>% rownames()
  # errorsMetastases <- predictionsMMAverage[metastases,] %>% filter(originalCall != predict) %>% nrow()
  # correctMetastases <- length(metastases) - errorsMetastases
  #
  # primary <- metaDataRef %>% filter(Status == "Primary") %>% rownames()
  # errorsPrimary <- predictionsMMAverage[primary, ] %>% filter(originalCall != predict) %>% nrow()
  # correctPrimary <- length(primary) - errorsPrimary

  #labelTypes <- c("Primary", "Recurrence", "Metastasis")

  for (i in seq(1:length(labelTypes))) {
    samples <-  metaDataRef %>% filter(!!sym(columnOfInterest) == labelTypes[i],
                                       rownames(.) %in% rownames(predictionsMM)) %>% rownames()
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

    #if(testSet == T) {
   #   samples <-  metaDataRef %>% filter(Status == labelTypes[i],
   #                                      rownames(.) %notin% throwOut) %>% rownames()
    #  samplesFiltered <-  predictionsMM[samples,] %>% filter(probability1 > probabilityScore) %>% rownames()
    #  correctSamplesThrowout <- predictionsMM[samplesFiltered,] %>% filter(originalCall == predict) %>% nrow()
    #  sampleTypeDF$errorsSamplesThrowout <- predictionsMM[samplesFiltered,] %>% filter(originalCall != predict) %>% nrow()
    #  sampleTypeDF$precisionThrowout <- round(correctSamplesThrowout / (correctSamplesThrowout + sampleTypeDF$errorsSamplesThrowout), digits = 3)
    #  sampleTypeDF$recallThrowout <- round(length(samplesFiltered) / length(samples), digits = 2)
    #}
    if (i == 1) {
      sampleTypeDFTotal <- sampleTypeDF

    } else {
      sampleTypeDFTotal <- rbind(sampleTypeDFTotal,
                                 sampleTypeDF)
    }

  }



  if (includeAllRow == T) {
    # Now add the results if all samples are taken into account
    #subsetResultTest <- predictionsMM %>% filter(originalCall %in% selectionMoreThanNCases)

    samples <-  metaDataRef %>% rownames()
    errorsSamples <- predictionsMM[samples,] %>% filter(originalCall != predict) %>% nrow()
    correctSamples <- length(samples) - errorsSamples
    samplesFiltered <- predictionsMM[samples,] %>% filter(probability1 > probabilityScore) %>% rownames()
    errorsFiltered <- predictionsMM[samplesFiltered,] %>% filter(originalCall != predict)%>% nrow()
    correctSamplesFiltered <- length(samplesFiltered) - errorsFiltered
    recall <- round(length(samplesFiltered) / length(samples), digits = 2)
  allSampleTypeDF <- data.frame(labelType = "All",
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


  sampleTypeDFTotal <- rbind(sampleTypeDFTotal, allSampleTypeDF)
}

 return(sampleTypeDFTotal)
}
