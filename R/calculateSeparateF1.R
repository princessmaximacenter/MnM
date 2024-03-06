#' Calculate performance measures for separate tumor (sub)types
#'
#' @param nSeeds How many seeds was the cross-validation setup run with?
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param crossValidation Specify whether the results are from the cross-validation setup or not.
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param nModels How many models were created for the majority voting system?
#' @param subtype  Do you want to obtain the predictions on the tumor subtype classification level?
#' @param throwOut  Are there samples you would like to remove from the test set due to poor data quality? If so, add their rownames here.
#' @param metaDataTest  Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the test set.
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the reference cohort.
#' @param filterOrNot Do you want to analyse the confident classifications only?
#' @param probabilityThreshold What is the threshold you would like to use to call a classification 'confident'?
#'
#' @return Dataframe containing the mean precision ($meanPrecision), F1-score ($meanF1), recall ($meanRecall)
#' and sensitivity per tumor (sub)type ($tumorType).
#' Results are stratified by the population frequency ($nCases).
#' @export
#'
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
    filterOrNot,
    probabilityThreshold) {




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
                                                    probabilityThreshold = probabilityThreshold,
                                                    filterOrNot = filterOrNot
                                                    )
    } else {
      fractionsCorrect <- extractIndividualValuesF1(predictionsMMFinal,
                                                    metaDataRef = metaDataRef,
                                                    classColumn = higherClassColumn,
                                                    probabilityThreshold = probabilityThreshold,
                                                    filterOrNot = filterOrNot)

    }
    fractionsCorrect$seed <- i
    if (i == 1) {
      accuracyDF <- fractionsCorrect
    } else {
      accuracyDF <- rbind(accuracyDF, fractionsCorrect)
    }
  }
  accuracyDF$nCases <- factor(accuracyDF$nCases, levels = unique(accuracyDF$nCases))
  meanNumbers <- accuracyDF %>% group_by(nCases, tumorType) %>%
    summarise(

      meanPrecision = mean(Precision),
      meanF1 = mean(F1),
      meanRecall = mean(Recall),
      meanSensitivity = mean(Sensitivity)
    )

 # levelsNCases <- c("n = 3", paste( nCases[-length(nCases)],"< n <=",nCases[-1])[-1], "n > 100", "All")
  #levelsNCases <- unique(meanNumbers$nCases)
  #meanNumbers$nCases <- factor(meanNumbers$nCases, levels = levelsNCases)
return(meanNumbers)
}
