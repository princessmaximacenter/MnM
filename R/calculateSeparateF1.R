#' Calculate performance measures for separate tumor (sub)types
#'
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param subtype  Do you want to obtain the classifications on the tumor subtype classification level?
#' @param metaDataTest  Metadata file containing the links between the samples and the tumor (sub)type diagnosis within the test set
#' @param probabilityThreshold What is the threshold you would like to use to call a classification 'confident'?
#' @param filterOrNot Do you want to filter the 'confident' classifications only for your calculation?
#'
#' @return Dataframe containing the mean precision ($meanPrecision), F1-score ($meanF1), recall ($meanRecall)
#' and sensitivity per tumor (sub)type ($tumorType).
#' Results are stratified by the population frequency ($nCases).
#' Lastly, it's specified whether a cross-validation ($type = Train) or test ($type = Test) type was used,
#' and whether the calculations were performed on the tumor type ($subtype = F) or subtype level ($subtype = T).
#' @export
#'
calculateSeparateF1 <- function(
    minorityDir,
    majorityDir,
    subtype,
    metaDataTest = NA,
    probabilityThreshold,
    filterOrNot = T) {



  `%notin%` <- base::Negate(`%in%`)
  allDirsMinority <- base::list.dirs(minorityDir, recursive = F)
  allDirsMajority <- base::list.dirs(majorityDir, recursive = F)
  selectedDirsMinority <- allDirsMinority[base::grep("seed", allDirsMinority)]
  selectedDirsMajority <- allDirsMajority[base::grep("seed", allDirsMajority)]

  if (base::length(selectedDirsMinority) != base::length(selectedDirsMajority)) {
    base::stop(base::paste("The number of models for the minority and majority classifier are not the same.",
         "Please check your models within the minorityDir and majorityDir"))
  } else if (!base::identical(base::sub(majorityDir, "", selectedDirsMajority), base::sub(minorityDir, "", selectedDirsMinority)) ) {
    base::stop(base::paste("It seems that classifications from the minority and majority classifier have not been run using the same seed.",
               "Please make sure you run the crossvalidation with the same seed for complementary classifications."))
  }

  if (base::length(base::grep("crossValidation", base::list.files(paste0(selectedDirsMajority, "/")))) > 0) {
    crossValidation <- T
  } else {
    crossValidation <- F
  }

  for (i in base::seq(1:base::length(selectedDirsMajority))) {
    if (crossValidation == T) {
    minorityDoc <- base::paste0(selectedDirsMinority[i],"/crossValidationMinorityResults.rds")
    majorityDoc <- base::paste0(selectedDirsMajority[i],"/crossValidationMajorityResults.rds")
    } else {
      minorityDoc <- base::paste0(minorityDir, "/minorityClassifierResult.rds")
      majorityDoc <- base::paste0(majorityDir, "/majorityClassifierResult.rds")

    }
    minority <- base::readRDS(minorityDoc)
    majority <- base::readRDS(majorityDoc)

    if (i == 1) {
    classColumn <- minority$metaDataRun$classColumn
    higherClassColumn <- minority$metaDataRun$higherClassColumn

    if (!base::is.na(metaDataTest)[1]) {
      base::print("Checking the performance for the test set based on values provided in dataframe 'metaDataTest'.")

      if (classColumn %notin% base::colnames(metaDataTest)) {
        base::print("Please note that the wanted column for the tumor subtype labels cannot be found within 'metaDataTest'.")
        base::print(base::paste0("Either change the column with the tumor subtype labels to the name: ", classColumn))
        base::stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
      } else if (higherClassColumn %notin% base::colnames(metaDataTest)) {

        base::print("Please note that the wanted column for the tumor type labels cannot be found within 'metaDataTest'.")
        base::print(base::paste0("Either change the column with the tumor type labels to the name: ", higherClassColumn))
        base::stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
      } else {
        base::print(base::paste0("Found columns ", classColumn, " and ", higherClassColumn, " within metaDataTest specifying the tumor subtype, and tumor type."))

        base::print("No original call found, adding it from metaDataTest")
      }
    }

    }
    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          subtype = subtype)

    predictionsMMFinal <- predictionsMMFinalList$predictionsMMFinal
    if (crossValidation == F ) {
      if (subtype == F) {
        predictionsMMFinal$originalCall <- metaDataTest[base::rownames(predictionsMMFinal), higherClassColumn]
      } else {
        predictionsMMFinal$originalCall <- metaDataTest[base::rownames(predictionsMMFinal), classColumn]
      }
    }

    if (subtype == T) {
      fractionsCorrect <- extractIndividualValuesF1(predictionsMM = predictionsMMFinal,
                                                    metaDataRef = minority$metaDataRef,
                                                    classColumn = classColumn,
                                                    probabilityThreshold = probabilityThreshold,
                                                    filterOrNot = filterOrNot
                                                    )
    } else {
      fractionsCorrect <- extractIndividualValuesF1(predictionsMM = predictionsMMFinal,
                                                    metaDataRef = minority$metaDataRef,
                                                    classColumn = higherClassColumn,
                                                    probabilityThreshold = probabilityThreshold,
                                                    filterOrNot = filterOrNot
                                                    )

    }
    fractionsCorrect$seed <- i
    if (i == 1) {
      accuracyDF <- fractionsCorrect
    } else {
      accuracyDF <- base::rbind(accuracyDF, fractionsCorrect)
    }
  }
  accuracyDF$nCases <- base::factor(accuracyDF$nCases, levels = base::unique(accuracyDF$nCases))
  meanNumbers <- accuracyDF %>%
    dplyr::group_by(nCases, tumorType) %>%
    dplyr::summarise(
      meanPrecision = base::mean(Precision),
      meanF1 = base::mean(F1),
      meanRecall = base::mean(Recall),
      meanSensitivity = base::mean(Sensitivity)
    )

  if (crossValidation == T) {
    meanNumbers$type <- "Train"
  } else {
    meanNumbers$type <- "Test"
  }

  meanNumbers$subtype <- subtype

 # levelsNCases <- c("n = 3", paste( nCases[-length(nCases)],"< n <=",nCases[-1])[-1], "n > 100", "All")
  #levelsNCases <- unique(meanNumbers$nCases)
  #meanNumbers$nCases <- factor(meanNumbers$nCases, levels = levelsNCases)
return(meanNumbers)
}
