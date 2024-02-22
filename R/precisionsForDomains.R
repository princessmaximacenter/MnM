precisionsForDomains <- function( nSeeds,
          classColumn,
          higherClassColumn,
          domainCol = "Domain",
          crossValidation,
          minorityDir,
          majorityDir,
          nModels,
          subtype,
          throwOut,
          metaDataTest,
          metaDataRef) {


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
      metaDataTest %<>% filter(rownames(.) %notin% throwOut)
      if (subtype == F) {
        predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), higherClassColumn]
      } else {
        predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), classColumn]
      }
    }





    if (subtype == T) {
      probabilityThreshold <- 0.7
    } else {
      probabilityThreshold <- 0.8

    }
    ourDomains <- unique(metaDataRef[,domainCol])

    if (crossValidation == F) {
      metaData <- metaDataTest
    } else {
      metaData <- metaDataRef
    }

    for (j in seq(1:length(ourDomains))) {

      samplesDomain <- metaData %>% filter(Domain == ourDomains[j]) %>% rownames(.)
      nSamples <- length(samplesDomain)
      predictionsMMDomain <- predictionsMMFinal %>% filter(rownames(.) %in% samplesDomain,
                                                           probability1 > probabilityThreshold)
      correct <- predictionsMMDomain %>% filter(originalCall == predict) %>% nrow()
      correctAll <- predictionsMMFinal %>% filter(rownames(.) %in% samplesDomain,
                                                  originalCall == predict) %>% nrow()
      accuracy <-  correctAll / length(samplesDomain)
      precision <- correct / nrow(predictionsMMDomain)
      recall <- nrow(predictionsMMDomain) / length(samplesDomain)

      ourDF <- data.frame(Domain = ourDomains[j],
                          nSamples = nSamples,
                          accuracy = accuracy,
                          precision = precision,
                          recall = recall)
      if (j == 1) {
        totalDF <- ourDF
      } else {
        totalDF <- rbind(totalDF, ourDF)

      }
    }

      nSamples <- nrow(predictionsMMFinal)
      predictionsMMFinalFiltered <- predictionsMMFinal %>% filter( probability1 > probabilityThreshold)
      correct <- predictionsMMFinalFiltered %>% filter(originalCall == predict) %>% nrow()
      correctAll <- predictionsMMFinal %>% filter(originalCall == predict) %>% nrow()
      accuracy <- correctAll / nrow(predictionsMMFinal)
      precision <- correct / nrow(predictionsMMFinalFiltered)
      recall <- nrow(predictionsMMFinalFiltered) / nrow(predictionsMMFinal)

      ourDF <- data.frame(Domain = "All",
                          nSamples = nSamples,
                          accuracy = accuracy,
                          precision = precision,
                          recall = recall)

      totalDF <- rbind(totalDF, ourDF)




      totalDF$seed <- i
    if (i == 1) {
      accuracyDF <- totalDF
    } else {
      accuracyDF <- rbind(accuracyDF, totalDF)
    }
  }

  meanNumbers <- accuracyDF %>% group_by(Domain) %>%
    summarise(
      nSamples = mean(nSamples),
      meanAccuracy = mean(accuracy),
      meanPrecision = mean(precision),
      meanRecall = mean(recall),
      sdAccuracy = sd(accuracy),
      sdPrecision = sd(precision),
      sdRecall = sd(recall)
    )

 return(meanNumbers)

}
