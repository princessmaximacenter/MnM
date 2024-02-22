calculateMeanAndSDOverall <- function(classColumn,
                                      higherClassColumn,
                                      minorityDir,
                                      majorityDir,
                                      metaDataRef,
                                      metaDataTest,
                                      throwOut = NA,
                                      subtype = F,
                                      crossValidation = T,
                                      nModels,
                                      nSeeds,
                                      subset = F,
                                      probabilityScoreThreshold

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
    if (crossValidation == F) {
      predictionsMMFinal %<>% filter(rownames(.) %notin% throwOut)
      if (subtype == F) {
        predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), higherClassColumn]
      } else {
        predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), classColumn]
      }
    }
    if (subset[1] != F) {
      predictionsMMFinal %<>% filter(rownames(.) %in% subset)
    }


      correct <- predictionsMMFinal %>% filter(originalCall == predict) %>% nrow()
      predictionsMMFiltered <- predictionsMMFinal %>% filter(probability1 > probabilityScoreThreshold)
      correctFiltered <- predictionsMMFiltered %>% filter(originalCall == predict) %>% nrow()
      accuracy <- correct / nrow(predictionsMMFinal)
      precision <- correctFiltered / nrow(predictionsMMFiltered)
      recall <- nrow(predictionsMMFiltered) / nrow(predictionsMMFinal)

      fractionsCorrect <- data.frame(accuracy = accuracy,
                                     precision = precision,
                                     recall = recall)

    fractionsCorrect$seed <- i
    if (i == 1) {
      accuracyDF <- fractionsCorrect
    } else {
      accuracyDF <- rbind(accuracyDF, fractionsCorrect)
    }
  }

  meanNumbers <- accuracyDF %>% #group_by(nCases) %>%
    summarise(
      meanAccuracy = mean(accuracy),
      meanPrecision = mean(precision),
      meanRecall = mean(recall),
      sdAccuracy = sd(accuracy),
      sdPrecision = sd(precision),
      sdRecall = sd(recall)
    )

  return(meanNumbers)

}
