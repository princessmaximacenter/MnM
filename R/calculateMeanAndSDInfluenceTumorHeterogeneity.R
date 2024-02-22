calculateMeanAndSDInfluenceTumorHeterogeneity <- function(
    minorityDir,
    majorityDir,
    trainOrTest,
    nSeeds,
    nModels,
    metaDataRef,
    metaDataTest = NA,
    classColumn,
    higherClassColumn,
    probabilityScoreTumor,
    probabilityScoreSubtype,
    throwOut,
    testSet
  ) {


  for (i in seq(1:nSeeds)) {
    if (trainOrTest == "Train") {
      crossValidation <- T
      minorityDoc <- paste0(minorityDir, "seed",i, "/crossValidationMinorityResults.rds")
      majorityDoc <- paste0(majorityDir, "seed",i, "/crossValidationMajorityResults.rds")
    } else {
      crossValidation <- F
      minorityDoc <- minorityDir
      majorityDoc <- majorityDir
    }


    minority <- readRDS(minorityDoc)
    majority <- readRDS(majorityDoc)

    if (trainOrTest == "Train") {
      crossValidation <- T
    } else {
      crossValidation <- F
    }
    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          nModels = nModels,
                                          subtype = F,
                                          metaDataRef = metaDataRef,
                                          classColumn = classColumn,
                                          higherClassColumn = higherClassColumn,
                                          crossValidation = crossValidation
    )

    predictionsMM <- predictionsMMFinalList$predictionsMMFinal

    #predictionsMM %<>% filter(rownames(.) %notin% throwOut)


    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          nModels = nModels,
                                          subtype = T,
                                          metaDataRef = metaDataRef,
                                          classColumn = classColumn,
                                          higherClassColumn = higherClassColumn,
                                          crossValidation = crossValidation
    )

    predictionsMMSubtype <- predictionsMMFinalList$predictionsMMFinal
    if (crossValidation == T) {
      predictionsMM$originalCall <- metaDataRef[rownames(predictionsMM), higherClassColumn]
    predictionsMMSubtype$originalCall <- metaDataRef[rownames(predictionsMMSubtype), classColumn]
    } else {
      predictionsMM$originalCall <- metaDataTest[rownames(predictionsMM), higherClassColumn]
      predictionsMMSubtype$originalCall <- metaDataTest[rownames(predictionsMMSubtype), classColumn]
    }

   if (trainOrTest == "Train") {
    statusDFLonger <- InfluenceTumorHeterogeneity(metaDataRef = metaDataRef,
                                                  predictionsMM = predictionsMM,
                                                  predictionsMMSubtype = predictionsMMSubtype,
                                                  trainOrTest = trainOrTest,
                                                  probabilityScoreTumor = probabilityScoreTumor,
                                                  probabilityScoreSubtype = probabilityScoreSubtype)
   } else if (trainOrTest == "Test") {
     predictionsMMSubtype %<>% filter(rownames(.) %notin% throwOut)
     predictionsMM %<>% filter(rownames(.) %notin% throwOut)
     statusDFLonger <- InfluenceTumorHeterogeneity(metaDataRef = metaDataTest,
                                                   predictionsMM = predictionsMM,
                                                   predictionsMMSubtype = predictionsMMSubtype,
                                                   trainOrTest = trainOrTest,
                                                   probabilityScoreTumor = probabilityScoreTumor,
                                                   probabilityScoreSubtype = probabilityScoreSubtype,
                                                   testSet = testSet)
   }
    statusDFLonger$seed <- i
    if (i == 1) {
      statusDFLongerTotal <- statusDFLonger
    } else {
      statusDFLongerTotal <- rbind(statusDFLongerTotal, statusDFLonger)
    }
  }
  statusDFLongerTotal
  meanNumbers <- statusDFLongerTotal %>% group_by(labelType,
                                                  trainOrTest,
                                                  type,
                                                  numberSamples
                                                  ) %>%
    summarise(
      meanAccuracy = mean(accuracy),
      meanPrecision = mean(precision),
      meanRecall = mean(recall),
      #meanIncorrect = mean(1 - fractionCorrect),
      #meanFractionIncorrectFiltered = mean(1 - fractionCorrectFiltered),
      sdAccuracy = sd(accuracy),
      sdPrecision = sd(precision),
      sdRecall = sd(recall)
    )

  return(meanNumbers)

}
