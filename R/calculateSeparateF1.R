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



  `%notin%` <- Negate(`%in%`)
  allDirsMinority <- list.dirs(minorityDir, recursive = F)
  allDirsMajority <- list.dirs(majorityDir, recursive = F)
  selectedDirsMinority <- allDirsMinority[grep("seed", allDirsMinority)]
  selectedDirsMajority <- allDirsMajority[grep("seed", allDirsMajority)]

  if (length(selectedDirsMinority) != length(selectedDirsMajority)) {
    stop("The number of models for the minority and majority classifier are not the same.
         Please check your models within the minorityDir and majorityDir")
  } else if (!identical(sub(majorityDir, "", selectedDirsMajority), sub(minorityDir, "", selectedDirsMinority)) ) {
    stop(paste("It seems that classifications from the minority and majority classifier have not been run using the same seed.",
               "\nPlease make sure you run the crossvalidation with the same seed for complementary classifications."))
  }

  if (length(grep("crossValidation", list.files(paste0(selectedDirsMajority, "/")))) > 0) {
    crossValidation <- T
  } else {
    crossValidation <- F
  }

  for (i in seq(1:length(selectedDirsMajority))) {
    if (crossValidation == T) {
    minorityDoc <- paste0(selectedDirsMinority[i],"/crossValidationMinorityResults.rds")
    majorityDoc <- paste0(selectedDirsMajority[i],"/crossValidationMajorityResults.rds")
    } else {
      minorityDoc <- paste0(minorityDir, "/minorityClassifierResult.rds")
      majorityDoc <- paste0(majorityDir, "/majorityClassifierResult.rds")

    }
    minority <- readRDS(minorityDoc)
    majority <- readRDS(majorityDoc)

    if (i == 1) {
    classColumn <- minority$metaDataRun$classColumn
    higherClassColumn <- minority$metaDataRun$higherClassColumn

    if (!is.na(metaDataTest)[1]) {
      print("Checking the performance for the test set based on values provided in dataframe 'metaDataTest'.")

      if (classColumn %notin% colnames(metaDataTest)) {
        print("Please note that the wanted column for the tumor subtype labels cannot be found within 'metaDataTest'.")
        print(paste0("Either change the column with the tumor subtype labels to the name: ", classColumn))
        stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
      } else if (higherClassColumn %notin% colnames(metaDataTest)) {

        print("Please note that the wanted column for the tumor type labels cannot be found within 'metaDataTest'.")
        print(paste0("Either change the column with the tumor type labels to the name: ", higherClassColumn))
        stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
      } else {
        print(paste0("Found columns ", classColumn, " and ", higherClassColumn, " within metaDataTest specifying the tumor subtype, and tumor type."))

        print("No original call found, adding it from metaDataTest")
      }
    }

    }
    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          subtype = subtype)

    predictionsMMFinal <- predictionsMMFinalList$predictionsMMFinal
    if (crossValidation == F ) {
      if (subtype == F) {
        predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), higherClassColumn]
      } else {
        predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), classColumn]
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
      accuracyDF <- rbind(accuracyDF, fractionsCorrect)
    }
  }
  accuracyDF$nCases <- factor(accuracyDF$nCases, levels = unique(accuracyDF$nCases))
  meanNumbers <- accuracyDF %>%
    dplyr::group_by(nCases, tumorType) %>%
    dplyr::summarise(
      meanPrecision = mean(Precision),
      meanF1 = mean(F1),
      meanRecall = mean(Recall),
      meanSensitivity = mean(Sensitivity)
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
