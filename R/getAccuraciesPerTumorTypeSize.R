#' Title
#'
#' @param predictionsMM
#' @param metaDataRef
#'
#' @return
#' @export
#'
#' @examples
getAccuraciesPerTumorTypeSize <- function(predictionsMM, probabilityThreshold = 0.8) {
  nCases <- c(1,3,5,10,20,40,100)
  fractionCorrect <- c()
  fractionCorrectFiltered <- c()
  nSamples <- c()
  nSamplesFiltered <- c()
  patientsPerTumor <- table(predictionsMM$originalCall)

  # Need to add the tumors that are not in the metadata, but are in the test set

  for ( i in 1:(length(nCases))){

    if ( i == length(nCases)){
      # Select tumors within a certain block (1-3, 4-5, 6-10, etc.)
      moreThanNCases <-  patientsPerTumor > nCases[i]
    } else {
      moreThanNCases <- (patientsPerTumor > nCases[i] &
                           patientsPerTumor <= nCases[i+1])
    }
    # Find their names

    selectionMoreThanNCases <- names(patientsPerTumor[moreThanNCases])

    subsetResultTest <- predictionsMM %>% filter(originalCall %in% selectionMoreThanNCases)

    filteredResultTest <- subsetResultTest %>% filter(probability1 > probabilityThreshold)

    nSamples[i] <- nrow(subsetResultTest)
    nSamplesFiltered[i] <- nrow(filteredResultTest)

    fractionCorrectFiltered[i] <- (filteredResultTest %>%
                                     filter(predict == originalCall) %>% nrow) / nrow(filteredResultTest)
    fractionCorrect[i] <- (subsetResultTest %>%
                             filter(predict == originalCall) %>% nrow) / nrow(subsetResultTest)
  }

  nCases <- c("n = 3", paste( nCases[-length(nCases)],"< n <=",nCases[-1])[-1], "n > 100")

  fractionsCorrect <- data.frame(nSamples = nSamples,
                                 nSamplesFiltered = nSamplesFiltered,
                                 fraction = round(nSamplesFiltered / nSamples, 2),
                                 fractionCorrect = round(fractionCorrect,2),
                                 fractionCorrectFiltered = round(fractionCorrectFiltered,2),
                                 fractionIncorrect = 1 - round(fractionCorrect,2),
                                 notClassified = 1 - round(nSamplesFiltered / nSamples, 2),
                                 nCases = nCases,
                                 numbersClassifiedFiltered = paste0("N = ", nSamplesFiltered),
                                 numbersClassifiedTotal = paste0("N = ", nSamples))


  fractionsCorrect$nCases <- factor(fractionsCorrect$nCases, levels = unique(fractionsCorrect$nCases))

  fractionsCorrect$fractionClassifiedCorrect <- round(fractionsCorrect$fraction * fractionsCorrect$fractionCorrectFiltered, digits = 2)
  fractionsCorrect$fractionIncorrect <- 1 - fractionsCorrect$fractionClassifiedCorrect - fractionsCorrect$notClassified
  return(fractionsCorrect)
}
