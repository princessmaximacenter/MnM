calculateMeanAndSDAccuracy <- function(classColumn,
                                       higherClassColumn,
                                       minorityDir,
                                       majorityDir,
                                       metaDataRef,
                                       metaDataTest = NA,
                                       throwOut = NA,
         subtype = F,
         crossValidation = T,
         nModels,
         nSeeds,
         rounding,
         probabilityThreshold
      #   subset = F
         ) {

  for (i in seq(1:nSeeds)) {
    if (crossValidation == T & nSeeds > 1) {
      minorityDoc <- paste0(minorityDir, "seed",i, "/crossValidationMinorityResults.rds")
      majorityDoc <- paste0(majorityDir, "seed",i, "/crossValidationMajorityResults.rds")
    } else {
      minorityDoc <- minorityDir
      majorityDoc <- majorityDir
    }

    minority <- readRDS(minorityDoc)
    majority <- readRDS(majorityDoc)

    predictionsMMFinalList <- integrateMM(minority = minority,
                                      majority = majority,
                                      nModels = nModels,
                                      subtype = subtype,
                                      metaDataRef = metaDataRef,
                                      classColumn = classColumn,
                                      higherClassColumn = higherClassColumn,
                                      crossValidation = crossValidation
    )

    predictionsMMFinal <- predictionsMMFinalList$predictionsMMFinal
    if (crossValidation == F ) {
      predictionsMMFinal %<>% filter(rownames(.) %notin% throwOut)
      if (subtype == F) {
      predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), higherClassColumn]
      } else {
      predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), classColumn]
    }
    }
   # if (subset[1] != F) {
   #   predictionsMMFinal %<>% filter(rownames(.) %in% subset)


   # }




    if (subtype == T) {
      fractionsCorrect <- getAccuraciesPerTumorTypeSize(predictionsMMFinal,
                                                        metaDataRef = metaDataRef,
                                                        classColumn = classColumn,
                                                        rounding = rounding,
                                                        probabilityThreshold = probabilityThreshold)
    } else {
      fractionsCorrect <- getAccuraciesPerTumorTypeSize(predictionsMMFinal,
                                                        metaDataRef = metaDataRef,
                                                        classColumn = higherClassColumn,
                                                        rounding = rounding,
                                                        probabilityThreshold = probabilityThreshold)

    }
    fractionsCorrect$seed <- i
    if (i == 1) {
      accuracyDF <- fractionsCorrect
    } else {
      accuracyDF <- rbind(accuracyDF, fractionsCorrect)
    }
  }

  meanNumbers <- accuracyDF %>% group_by(nCases) %>%
    summarise(
      meanFractionCorrect = mean(fractionCorrect),
      meanFractionCorrectFiltered = mean(fractionCorrectFiltered),
      meanFractionIncorrect = mean(1 - fractionCorrect),
      meanFractionIncorrectFiltered = mean(1 - fractionCorrectFiltered),
      sdFractionCorrect = sd(fractionCorrect),
      sdFractionCorrectFiltered = sd(fractionCorrectFiltered),
      meanCasesFiltered = round(mean(nSamplesFiltered), digits = 0),
      meanPrecision = mean(Precision),
      meanF1 = mean(F1),
      meanFractionCorrect2 = mean(fractionCorrect2),
      meanFractionCorrect3 = mean(fractionCorrect3),
      sdFractionCorrect2 = sd(fractionCorrect2),
      sdFractionCorrect3 = sd(fractionCorrect3),
      meanRecall = mean(Recall),
      medianF1 = median(F1),
      sdPrecision = sd(Precision),
      sdF1 = sd(F1),
      sdRecall = sd(Recall)
    )

  return(meanNumbers)
}
