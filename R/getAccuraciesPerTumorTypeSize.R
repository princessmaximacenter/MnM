#' Obtain accuracies per tumor (sub)type sample size
#'
#' @param predictionsMM Dataframe showing the top 3 predictions for the tumor (sub)type, together with their probability scores.
#' @param metaDataRef Metadata-file for the reference cohort containing the tumor (sub)type labels.
#' @param classColumn Which column within the metadata file contains the tumor (sub)type labels?
#' @return Dataframe showing the accuracies for the different tumor size
#' @export
#'
getAccuraciesPerTumorTypeSize <- function(predictionsMM,
                                          metaDataRef,
                                          classColumn,
                                          rounding = T,
                                          probabilityThreshold) {
  nCases <- c(1,5,10,20,40,100)
  fractionCorrect <- c()
  fractionCorrect2 <- c()
  fractionCorrect3 <- c()
  fractionCorrectFiltered <- c()
  nSamples <- c()
  nSamplesFiltered <- c()
  errorsFiltered <- c()
  patientsPerTumor <- table(metaDataRef[,classColumn])
  sensitivity <- c()
  specificity <- c()
  Precision <- c()
  Recall <- c()
  F1 <- c()

  predictionsMMFiltered <- predictionsMM %>% filter(probability1 > probabilityThreshold)
  tumorConfusionMatrix <- confusionMatrix(factor(predictionsMMFiltered$predict,
                                                 levels = unique(c(predictionsMMFiltered$originalCall))),
                                          factor(predictionsMMFiltered$originalCall, levels =  unique(c(predictionsMMFiltered$originalCall))),
                                          dnn = c("Prediction", "Reference"))
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

    tumorConfusionMatrixSelection <- tumorConfusionMatrix$byClass %>% as.data.frame() %>% filter(rownames(.) %in% paste0("Class: ", selectionMoreThanNCases))

    tumorConfusionMatrixSelection$F1[is.na(tumorConfusionMatrixSelection$F1)] <- 0
    tumorConfusionMatrixSelection$Precision[is.na(tumorConfusionMatrixSelection$Precision)] <- 0
    recall <- table(filteredResultTest[,"originalCall"]) / patientsPerTumor[names(patientsPerTumor) %in% names(table(filteredResultTest[,"originalCall"]))]

    nSamples[i] <- nrow(subsetResultTest)
    nSamplesFiltered[i] <- nrow(filteredResultTest)

    fractionCorrectFiltered[i] <- (filteredResultTest %>%
                                     filter(predict == originalCall) %>% nrow) / nrow(filteredResultTest)
    errorsFiltered[i] <- filteredResultTest %>%
                            filter(predict != originalCall) %>% nrow
    fractionCorrect[i] <- (subsetResultTest %>%
                             filter(predict == originalCall) %>% nrow) / nrow(subsetResultTest)
    fractionCorrect2[i] <- (subsetResultTest %>%
                           filter(predict == originalCall | predict2 == originalCall) %>% nrow) / nrow(subsetResultTest)
    fractionCorrect3[i] <- (subsetResultTest %>%
                           filter(predict == originalCall | predict2 == originalCall |
                                    predict3 == originalCall) %>% nrow) / nrow(subsetResultTest)


    #sensitivity[i] <- mean(tumorConfusionMatrixSelection$sensitivity)
    #specificity[i] <- mean(tumorConfusionMatrixSelection$specificity)
    Precision[i] <- mean(tumorConfusionMatrixSelection$Precision)
    #Recall[i] <- mean(tumorConfusionMatrixSelection$Recall)
    F1[i] <- mean(tumorConfusionMatrixSelection$F1)
    Recall[i] <- mean(recall)# Make this how many of the samples are included into the classified samples

  }

  #nCases <- c("n = 3", paste0("n = ", nCases[-length(nCases)] + 1,"-",nCases[-1])[-1], "n > 100")
  #nCases <- c("n = 3", paste( nCases[-c(1,length(nCases))],"< n <=",nCases[-c(1,2)])[-1], "n > 100")
  nCases  <- c("3 < n <= 5", paste( nCases[-c(length(nCases))],"< n <=",nCases[-c(1)])[-1], "n > 100")
  if (rounding == T) {
    fractionsCorrect <- data.frame(nSamples = nSamples,
                                                     nSamplesFiltered = nSamplesFiltered,
                                                     fraction = round(nSamplesFiltered / nSamples, 2),
                                                     fractionCorrect = round(fractionCorrect,2),
                                                     fractionCorrect2 = round(fractionCorrect2, 2),
                                                     fractionCorrect3 = round(fractionCorrect3, 2),
                                                     fractionCorrectFiltered = round(fractionCorrectFiltered,2),
                                                     fractionIncorrect = 1 - round(fractionCorrect,2),
                                                     notClassified = 1 - round(nSamplesFiltered / nSamples, 2),
                                                     nCases = nCases,
                                                     numbersClassifiedFiltered = paste0("N = ", nSamplesFiltered),
                                                     numbersClassifiedTotal = paste0("N = ", nSamples),
                                                     errorsFiltered = paste0("N = ", errorsFiltered),
                                                     Precision = Precision,
                                   F1 = F1,
                                   Recall = Recall
                                   )
  } else {
    fractionsCorrect <- data.frame(nSamples = nSamples,
                                   nSamplesFiltered = nSamplesFiltered,
                                   fraction = nSamplesFiltered / nSamples,
                                   fractionCorrect = fractionCorrect,
                                   fractionCorrect2 =fractionCorrect2,
                                   fractionCorrect3 = fractionCorrect3,
                                   fractionCorrectFiltered = fractionCorrectFiltered,
                                   fractionIncorrect = 1 - fractionCorrect,
                                   notClassified = 1 - nSamplesFiltered / nSamples,
                                   nCases = nCases,
                                   numbersClassifiedFiltered = paste0("N = ", nSamplesFiltered),
                                   numbersClassifiedTotal = paste0("N = ", nSamples),
                                   errorsFiltered = paste0("N = ", errorsFiltered),
                                   Precision = Precision,
                                   F1 = F1,
                                   Recall = Recall)
    }

  # Now add the results if all samples are taken into account
  #subsetResultTest <- predictionsMM %>% filter(originalCall %in% selectionMoreThanNCases)

  filteredResultTest <- predictionsMM %>% filter(probability1 > probabilityThreshold)

  nSamples <- nrow(predictionsMM)
  nSamplesFiltered <- nrow(filteredResultTest)

  fractionCorrectFiltered <- (filteredResultTest %>%
                                   filter(predict == originalCall) %>% nrow) / nrow(filteredResultTest)
  errorsFiltered <- filteredResultTest %>%
    filter(predict != originalCall) %>% nrow
  fractionCorrect <- (predictionsMM %>%
                           filter(predict == originalCall) %>% nrow) / nrow(predictionsMM)
  fractionCorrect2 <- (predictionsMM %>%
                          filter(predict == originalCall | predict2 == originalCall) %>% nrow) / nrow(predictionsMM)
  fractionCorrect3 <- (predictionsMM %>%
                         filter(predict == originalCall | predict2 == originalCall |
                                predict3 == originalCall) %>% nrow) / nrow(predictionsMM)

  recall <- table(filteredResultTest[,"originalCall"]) / patientsPerTumor[names(patientsPerTumor) %in% names(table(filteredResultTest[,"originalCall"]))]
  tumorConfusionMatrixSelection <- tumorConfusionMatrix$byClass %>% as.data.frame()
  tumorConfusionMatrixSelection$Precision[is.na(tumorConfusionMatrixSelection$Precision)] <- 0
  tumorConfusionMatrixSelection$F1[is.na(tumorConfusionMatrixSelection$F1)] <- 0
  if (rounding == T) {
    allSamplesCorrect <-  data.frame(nSamples = nSamples,
                                     nSamplesFiltered = nSamplesFiltered,
                                     fraction = round(nSamplesFiltered / nSamples, 2),
                                     fractionCorrect = round(fractionCorrect,2),
                                     fractionCorrect2 = round(fractionCorrect2, 2),
                                     fractionCorrect3 = round(fractionCorrect3, 2),
                                     fractionCorrectFiltered = round(fractionCorrectFiltered,2),
                                     fractionIncorrect = 1 - round(fractionCorrect,2),
                                     notClassified = 1 - round(nSamplesFiltered / nSamples, 2),
                                     nCases = "All",
                                     numbersClassifiedFiltered = paste0("N = ", nSamplesFiltered),
                                     numbersClassifiedTotal = paste0("N = ", nSamples),
                                     errorsFiltered = paste0("N = ", errorsFiltered),
                                     Precision = mean(tumorConfusionMatrixSelection$Precision),
                                     F1= mean(tumorConfusionMatrixSelection$F1),
                                     Recall = mean(recall)
                                     )
  } else {
    allSamplesCorrect <-  data.frame(nSamples = nSamples,
                                      nSamplesFiltered = nSamplesFiltered,
                                      fraction = nSamplesFiltered / nSamples,
                                      fractionCorrect = fractionCorrect,
                                     fractionCorrect2 = fractionCorrect2,
                                     fractionCorrect3 = fractionCorrect3,
                                      fractionCorrectFiltered = fractionCorrectFiltered,
                                      fractionIncorrect = 1 - fractionCorrect,
                                      notClassified = 1 - nSamplesFiltered / nSamples,
                                      nCases = "All",
                                      numbersClassifiedFiltered = paste0("N = ", nSamplesFiltered),
                                      numbersClassifiedTotal = paste0("N = ", nSamples),
                                      errorsFiltered = paste0("N = ", errorsFiltered),
                                     Precision = mean(tumorConfusionMatrixSelection$Precision),
                                     F1 = mean(tumorConfusionMatrixSelection$F1),
                                     Recall = mean(recall))
  }


  fractionsCorrect <- rbind(fractionsCorrect, allSamplesCorrect)

  fractionsCorrect$nCases <- factor(fractionsCorrect$nCases, levels = unique(fractionsCorrect$nCases))
  if (rounding == T) {
  fractionsCorrect$fractionClassifiedCorrect <- round(fractionsCorrect$fraction * fractionsCorrect$fractionCorrectFiltered, digits = 2)
  } else {
    fractionsCorrect$fractionClassifiedCorrect <- fractionsCorrect$fraction * fractionsCorrect$fractionCorrectFiltered
  }
  fractionsCorrect$fractionClassifiedIncorrect <- 1 - fractionsCorrect$fractionClassifiedCorrect - fractionsCorrect$notClassified
  return(fractionsCorrect)
}
