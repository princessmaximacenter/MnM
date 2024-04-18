#' Compare results between samples from different domains
#'
#' @param nSeeds How many seeds was the cross-validation setup run with?
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param domainColumn Column in the metadata file that contains the domain labels.
#' @param crossValidation Specify whether the results are from the cross-validation setup or not.
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param nModels How many models were created for the majority voting system?
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' @param throwOut Are there samples you would like to remove from the test set due to poor data quality? If so, add their rownames here.
#' @param metaDataTest Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the test set.
#' @param probabilityThreshold What is the probability score threshold you would like to use to call a classification 'confident' for M&M?
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#'
#' @return Dataframe containing the mean accuracy ($meanAccuracy),
#' precision ($meanPrecision), and recall ($meanRecall),
#' together with their standard deviations ($sdAccuracy,
#' $sdPrecision, $sdRecall), and the total number of samples associated to the domain group ($nSamples)
#' @export
#'
precisionsForDomains <- function( nSeeds,
          classColumn,
          higherClassColumn,
          domainColumn,
          crossValidation,
          minorityDir,
          majorityDir,
          nModels,
          subtype,
          throwOut,
          metaDataTest,
          probabilityThreshold,
          metaDataRef) {

  `%notin%` <- Negate(`%in%`)
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
                                          subtype = subtype,
                                          classColumn = classColumn,
                                          higherClassColumn = higherClassColumn
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

    ourDomains <- unique(metaDataRef[,domainColumn])

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
