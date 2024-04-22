#' Get performance measures for separate heterogeneous groups within data set
#'
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the reference cohort.
#' @param predictionsMM Dataframe containing the top 3 final tumor type classification labels ($predict{2,3}) with their accompanying probability scores
#' ($probability{1,2,3} and the original diagnosis label ($originalCall).
#' @param probabilityScore What is the probability score threshold you would like to use to call a classification 'confident'?
#' @param columnOfInterest Which column contains the different groups that should be analyzed?
#' @param labelTypes What are the names for the different groups within the column of interest (e.g. Primary, Recurrence, Metastasis)?
#' @param includeAllRow Do you want to also generate an entry that specifies the performance for all samples taken together?
#'
#' @return Dataframe containing the total number of samples within the group ($numberSamples), how many of those were correct ($correctSamples),
#' how many errors were present ($errorsSamples),
#' what was the accuracy, how many samples were confidently classified ($numberSamplesFiltered),
#' how many of those confident classifications were correct ($errorsFiltered),
#' what was the accuracy within the confident classifications ($precision) and how many samples reached the confidence treshold ($recall).

primaryRecurrenceMetastasisNumbers <- function(metaDataRef,
         predictionsMM,
         probabilityScore,
         columnOfInterest,
         labelTypes,
         includeAllRow = F
         ) {


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

    if (i == 1) {
      sampleTypeDFTotal <- sampleTypeDF

    } else {
      sampleTypeDFTotal <- rbind(sampleTypeDFTotal,
                                 sampleTypeDF)
    }

  }



  if (includeAllRow == T) {
    # Now add the results if all samples are taken into account

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
