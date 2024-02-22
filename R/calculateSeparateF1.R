calculateSeparateF1 <- function(
    nSeeds,
    classColumn,
    higherClassColumn,
    crossValidation,
    minorityDir,
    majorityDir,
    nModels,
    subtype,
    throwOut,
    metaDataTest,
    metaDataRef,
    filterOrNot) {




  for (i in seq(1:nSeeds)) {
    if (crossValidation == T) {
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





    if (subtype == T) {
      fractionsCorrect <- extractIndividualValuesF1(predictionsMMFinal,
                                                    metaDataRef = metaDataRef,
                                                    classColumn = classColumn,
                                                    probabilityThreshold = 0.72,
                                                    filterOrNot = filterOrNot
                                                    )
    } else {
      fractionsCorrect <- extractIndividualValuesF1(predictionsMMFinal,
                                                    metaDataRef = metaDataRef,
                                                    classColumn = higherClassColumn,
                                                    probabilityThreshold = 0.8,
                                                    filterOrNot = filterOrNot)

    }
    fractionsCorrect$seed <- i
    if (i == 1) {
      accuracyDF <- fractionsCorrect
    } else {
      accuracyDF <- rbind(accuracyDF, fractionsCorrect)
    }
  }

  meanNumbers <- accuracyDF %>% group_by(tumorType, nCases) %>%
    summarise(

      meanPrecision = mean(Precision),
      meanF1 = mean(F1),
      meanRecall = mean(Recall),
      meanSensitivity = mean(Sensitivity)
    )

  nCases <- c(1,3,5,10,20,40,100)
 # levelsNCases <- c("n = 3", paste( nCases[-length(nCases)],"< n <=",nCases[-1])[-1], "n > 100", "All")
  levelsNCases <- c("n = 3", paste( nCases[-length(nCases)],"< n <=",nCases[-1])[-1], "n > 100")
  meanNumbers$nCases <- factor(meanNumbers$nCases, levels = levelsNCases)
return(meanNumbers)
}
