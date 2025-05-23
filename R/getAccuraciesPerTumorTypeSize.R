#' Obtain accuracy per tumor (sub)type sample size
#'
#' @param predictionsMM Dataframe showing the top 3 classifications for the tumor (sub)type, together with their probability scores.
#' @param metaDataRef Metadata-file for the reference cohort containing the tumor (sub)type labels.
#' @param classColumn Which column within the metadata file contains the tumor (sub)type labels?
#' @param rounding Do you want rounded numbers for the performance scores? Default is TRUE.
#' @param probabilityThreshold What is the probability score threshold you would like to use to call a classification 'confident'?
#' @param withoutF1 Do you need to run the calculations without the F1-score and precision?
#' This is the case when there are only 1-2 sample types within the dataset, caret::confusionMatrix crashes.
#' This parameter prevents the crashing.
#' @return Dataframe showing the performance scores for the different tumor frequencies (nCases).
#' Shown are the total amount of samples within the range ($nSamples),
#' how many are classified confidently ($nSamplesFiltered),
#' what is the fraction classified confidently ($fraction),
#' what is the fraction of correct and incorrect classifications within the highest scoring classification ($fractionCorrect and $fractionIncorrect),
#' the top 2  ($fractionCorrect2) an top 3 highest scoring classifications ($fractioncorrect3).
#' Furthermore, the fraction correct classifications within the confident classifications is denoted ($fractionCorrectFiltered),
#' together with the amount of errors remaining in here ($errorsFiltered).
#'
#' For the rest, for the confident classifications the average precision ($Precision),
#' F1 score ($F1) and recall ($Recall) for all tumor entities combined is determined.
#'@import caret
#'@export
getAccuraciesPerTumorTypeSize <- function(predictionsMM,
                                          metaDataRef,
                                          classColumn,
                                          rounding = T,
                                          probabilityThreshold,
                                          withoutF1= F) {
  nCases <- c(0, 5,10,20,40,100)
  fractionCorrect <- c()
  fractionCorrect2 <- c()
  fractionCorrect3 <- c()
  fractionCorrectFiltered <- c()
  nSamples <- c()
  nSamplesFiltered <- c()
  errorsFiltered <- c()
  patientsPerTumor <- base::table(metaDataRef[,classColumn])
  sensitivity <- c()
  specificity <- c()
  Precision <- c()
  Recall <- c()
  F1 <- c()
  predictionsMM <- predictionsMM %>% dplyr::filter(originalCall %in% base::unique(metaDataRef[,classColumn]))
  predictionsMMFiltered <- predictionsMM %>% dplyr::filter(probability1 > probabilityThreshold)

  if (withoutF1 == F) {
    tumorConfusionMatrix <- caret::confusionMatrix(base::factor(predictionsMMFiltered$predict,
                                                                levels = base::unique(c(predictionsMMFiltered$originalCall))),
                                                   factor(predictionsMMFiltered$originalCall, levels =  base::unique(c(predictionsMMFiltered$originalCall))),
                                                   dnn = c("Prediction", "Reference"))
  }

  # Need to add the tumors that are not in the metadata, but are in the test set

  for ( i in 1:(base::length(nCases))){

    if ( i == base::length(nCases)){
      # Select tumors within a certain block (1-3, 4-5, 6-10, etc.)
      moreThanNCases <-  patientsPerTumor > nCases[i]
    } else {
      moreThanNCases <- (patientsPerTumor > nCases[i] &
                           patientsPerTumor <= nCases[i+1])
    }
    # Find their names

    selectionMoreThanNCases <- base::names(patientsPerTumor[moreThanNCases])

    subsetResultTest <- predictionsMM %>% dplyr::filter(originalCall %in% selectionMoreThanNCases)

    filteredResultTest <- subsetResultTest %>% dplyr::filter(probability1 > probabilityThreshold)

    if (withoutF1 == F) {
      tumorConfusionMatrixSelection <- tumorConfusionMatrix$byClass %>% base::as.data.frame() %>%
        dplyr::filter(base::rownames(.) %in% base::paste0("Class: ", selectionMoreThanNCases))

      tumorConfusionMatrixSelection$F1[base::is.na(tumorConfusionMatrixSelection$F1)] <- 0
    }



    #tumorConfusionMatrixSelection$Precision[is.na(tumorConfusionMatrixSelection$Precision)] <- 0
    #recall <- base::table(filteredResultTest[,"originalCall"]) / patientsPerTumor[base::names(patientsPerTumor) %in% base::names(base::table(filteredResultTest[,"originalCall"]))]

    subsetResultTest$originalCall <- factor(subsetResultTest$originalCall)
    filteredResultTest$originalCall <- factor(filteredResultTest$originalCall, levels = levels(subsetResultTest$originalCall))

    #tumorConfusionMatrixSelection$Precision[is.na(tumorConfusionMatrixSelection$Precision)] <- 0
    recall <- base::table(filteredResultTest[,"originalCall"]) / base::table(subsetResultTest[, "originalCall"])

    #recall <- base::table(filteredResultTest[,"originalCall"]) / base::table(subsetResultTest[base::names(base::table(subsetResultTest[, "originalCall"])) %in% base::names(base::table(filteredResultTest[,"originalCall"])), "originalCall"])
    nSamples[i] <- base::nrow(subsetResultTest)
    nSamplesFiltered[i] <- base::nrow(filteredResultTest)

    fractionCorrectFiltered[i] <- (filteredResultTest %>%
                                     dplyr::filter(predict == originalCall) %>% base::nrow()) / base::nrow(filteredResultTest)
    errorsFiltered[i] <- filteredResultTest %>%
      dplyr::filter(predict != originalCall) %>% base::nrow()
    fractionCorrect[i] <- (subsetResultTest %>%
                             dplyr::filter(predict == originalCall) %>% base::nrow()) / base::nrow(subsetResultTest)
    fractionCorrect2[i] <- (subsetResultTest %>%
                              dplyr::filter(predict == originalCall | predict2 == originalCall) %>% base::nrow()) / nrow(subsetResultTest)
    fractionCorrect3[i] <- (subsetResultTest %>%
                              dplyr::filter(predict == originalCall | predict2 == originalCall |
                                              predict3 == originalCall) %>% base::nrow()) / base::nrow(subsetResultTest)

    if (withoutF1 == F) {
      Precision[i] <- base::mean(tumorConfusionMatrixSelection$Precision, na.rm = T)
      #Recall[i] <- mean(tumorConfusionMatrixSelection$Recall)
      F1[i] <- base::mean(tumorConfusionMatrixSelection$F1)
    }

    Recall[i] <- base::mean(recall)# Make this how many of the samples are included into the classified samples

  }

  nCases  <- c("3 <= n <= 5", base::paste( nCases[-c(base::length(nCases))],"< n <=",nCases[-c(1)])[-1], "n > 100")
  if (rounding == T) {
    fractionsCorrect <- data.frame(nSamples = nSamples,
                                   nSamplesFiltered = nSamplesFiltered,
                                   fraction = base::round(nSamplesFiltered / nSamples, 2),
                                   fractionCorrect = base::round(fractionCorrect,2),
                                   fractionCorrect2 = base::round(fractionCorrect2, 2),
                                   fractionCorrect3 = base::round(fractionCorrect3, 2),
                                   fractionCorrectFiltered = base::round(fractionCorrectFiltered,2),
                                   fractionIncorrect = 1 - base::round(fractionCorrect,2),
                                   notClassified = 1 - base::round(nSamplesFiltered / nSamples, 2),
                                   nCases = nCases,
                                   numbersClassifiedFiltered = base::paste0("N = ", nSamplesFiltered),
                                   numbersClassifiedTotal = base::paste0("N = ", nSamples),
                                   errorsFiltered = base::paste0("N = ", errorsFiltered),
                                   Recall = Recall
    )

    if (withoutF1 == F) {
      fractionsCorrect$Precision <- Precision
      fractionsCorrect$F1 <- F1
    }

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
                                   numbersClassifiedFiltered = base::paste0("N = ", nSamplesFiltered),
                                   numbersClassifiedTotal = base::paste0("N = ", nSamples),
                                   errorsFiltered = base::paste0("N = ", errorsFiltered),
                                   Recall = Recall)

    if (withoutF1 == F) {
      fractionsCorrect$Precision <- Precision
      fractionsCorrect$F1 <- F1
    }
  }

  # Now add the results if all samples are taken into account
  #subsetResultTest <- predictionsMM %>% filter(originalCall %in% selectionMoreThanNCases)

  filteredResultTest <- predictionsMM %>% dplyr::filter(probability1 > probabilityThreshold)

  nSamples <- base::nrow(predictionsMM)
  nSamplesFiltered <- base::nrow(filteredResultTest)

  fractionCorrectFiltered <- (filteredResultTest %>%
                                dplyr::filter(predict == originalCall) %>% base::nrow()) / base::nrow(filteredResultTest)
  errorsFiltered <- filteredResultTest %>%
    dplyr::filter(predict != originalCall) %>% base::nrow()
  fractionCorrect <- (predictionsMM %>%
                        dplyr::filter(predict == originalCall) %>% base::nrow()) / base::nrow(predictionsMM)
  fractionCorrect2 <- (predictionsMM %>%
                         dplyr::filter(predict == originalCall | predict2 == originalCall) %>% base::nrow() ) / base::nrow(predictionsMM)
  fractionCorrect3 <- (predictionsMM %>%
                         dplyr::filter(predict == originalCall | predict2 == originalCall |
                                         predict3 == originalCall) %>% base::nrow(.)) / base::nrow(predictionsMM)

  predictionsMM$originalCall <- factor(predictionsMM$originalCall)
  filteredResultTest$originalCall <- factor(filteredResultTest$originalCall, levels = levels(predictionsMM$originalCall))
  #recall <- base::table(filteredResultTest[,"originalCall"]) / patientsPerTumor[base::names(patientsPerTumor) %in% base::names(base::table(filteredResultTest[,"originalCall"]))]

  recall <- base::table(filteredResultTest[,"originalCall"]) / base::table(predictionsMM[, "originalCall"])

  if (withoutF1 == F) {
    tumorConfusionMatrixSelection <- tumorConfusionMatrix$byClass %>% base::as.data.frame()
    tumorConfusionMatrixSelection$F1[is.na(tumorConfusionMatrixSelection$F1)] <- 0
  }

  if (rounding == T) {
    allSamplesCorrect <-  data.frame(nSamples = nSamples,
                                     nSamplesFiltered = nSamplesFiltered,
                                     fraction = base::round(nSamplesFiltered / nSamples, 2),
                                     fractionCorrect = base::round(fractionCorrect,2),
                                     fractionCorrect2 = base::round(fractionCorrect2, 2),
                                     fractionCorrect3 = base::round(fractionCorrect3, 2),
                                     fractionCorrectFiltered = base::round(fractionCorrectFiltered,2),
                                     fractionIncorrect = 1 - base::round(fractionCorrect,2),
                                     notClassified = 1 - base::round(nSamplesFiltered / nSamples, 2),
                                     nCases = "All",
                                     numbersClassifiedFiltered = base::paste0("N = ", nSamplesFiltered),
                                     numbersClassifiedTotal = base::paste0("N = ", nSamples),
                                     errorsFiltered = base::paste0("N = ", errorsFiltered),
                                     Recall = base::mean(recall)
    )

    if (withoutF1 == F) {

      allSamplesCorrect$Precision <- base::mean(tumorConfusionMatrixSelection$Precision, na.rm = T)
      allSamplesCorrect$F1 <- base::mean(tumorConfusionMatrixSelection$F1)
    }

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
                                     numbersClassifiedFiltered = base::paste0("N = ", nSamplesFiltered),
                                     numbersClassifiedTotal = base::paste0("N = ", nSamples),
                                     errorsFiltered = base::paste0("N = ", errorsFiltered),
                                     Recall = base::mean(recall))

    if (withoutF1 == F) {

      allSamplesCorrect$Precision <- base::mean(tumorConfusionMatrixSelection$Precision, na.rm = T)
      allSamplesCorrect$F1 <- base::mean(tumorConfusionMatrixSelection$F1)
    }
  }


  fractionsCorrect <- base::rbind(fractionsCorrect, allSamplesCorrect)

  fractionsCorrect$nCases <- base::factor(fractionsCorrect$nCases, levels = base::unique(fractionsCorrect$nCases))
  if (rounding == T) {
    fractionsCorrect$fractionClassifiedCorrect <- base::round(fractionsCorrect$fraction * fractionsCorrect$fractionCorrectFiltered, digits = 2)
  } else {
    fractionsCorrect$fractionClassifiedCorrect <- fractionsCorrect$fraction * fractionsCorrect$fractionCorrectFiltered
  }
  fractionsCorrect$fractionClassifiedIncorrect <- 1 - fractionsCorrect$fractionClassifiedCorrect - fractionsCorrect$notClassified
  return(fractionsCorrect)
}
