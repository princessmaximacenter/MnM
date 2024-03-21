#' Calculate mean and standard deviation for whole dataset
#'
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the reference cohort.
#' @param metaDataTest Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the test set.
#' @param throwOut Are there samples you would like to remove from the test set due to poor data quality? If so, add their rownames here.
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' @param crossValidation Specify whether the results are from the cross-validation setup or not.
#' @param nModels How many models were created for the majority voting system?
#' @param nSeeds How many seeds was the cross-validation setup run with?
#' @param subset Do you want to only include a subset of the data for the performance measure calculations? If so, add their rownames here.
#' @param probabilityThreshold What is the threshold you would like to use to call a classification 'confident'?
#'
#' @return Dataframe with performance for all samples within the dataset combined,
#' instead of with stratification by tumor (sub)type frequency within reference cohort.
#'
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
                                      probabilityThreshold
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
    predictionsMMFiltered <- predictionsMMFinal %>% filter(probability1 > probabilityThreshold)
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
