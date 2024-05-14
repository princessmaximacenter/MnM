#' Calculate performance of M&M on total train and test set
#'
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param metaDataTest Metadata file containing the links between the samples and the tumor type and subtype diagnoses within the test set.
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level (subtype = TRUE)?
#' @param rounding Do you want rounded numbers for the performance scores? Default = FALSE.
#' @param probabilityThreshold What is the probability score threshold you would like to use to call a classification 'confident'?
#'
#' @return Dataframe containing the average performance of tumor classifications within a certain frequency range (nCases).
#' Included are the averages for the percentage of correctly and incorrectly classified samples ($meanFractionCorrect and $meanFractionIncorrect),
#' correctly and incorrectly classified 'confident' samples ($meanFractionCorrectFiltered and $meanFractionIncorrectFiltered),
#' averaged precision values for the average of all tumor entities within the frequency range ($meanPrecision),
#' averaged recall values for the average of all tumor entities within the frequency range ($meanRecall),
#' and averaged F1 scores for the average of all tumor entities within the frequency range ($meanF1).
#' Furthermore, the standard deviations of all mentioned values are calculated in case multiple runs with different seeds have been performed ($sd...).
#' Please note that the precision, F1 and recall are calculated for the confident sample classifications only.
#'The total amount of samples within each frequency range ($meanSamples) is also specified.
#' Lastly, it's specified whether a cross-validation (Train) or test (Test) type was used.
#' @export
#' @import stats
#'
calculateMeanAndSDAccuracy <- function(
  minorityDir,
  majorityDir,
  metaDataTest = NA,
         subtype = F,
         rounding = F,
         probabilityThreshold
         ) {

  `%notin%` <- base::Negate(`%in%`)
  if ("minorityClassifierResult.rds" %notin% base::list.files(minorityDir)) {
    crossValidation <- T
  allDirsMinority <- base::list.dirs(minorityDir, recursive = F)
  allDirsMajority <- base::list.dirs(majorityDir, recursive = F)
  selectedDirsMinority <- allDirsMinority[base::grep("seed", allDirsMinority)]
  selectedDirsMajority <- allDirsMajority[base::grep("seed", allDirsMajority)]

  base::print(base::paste0("Found ",base::length(selectedDirsMajority), " directories with different cross-validation runs.",
  " Calculating average performance values for all combined."))
  if (base::length(selectedDirsMinority) != base::length(selectedDirsMajority)) {
    base::stop("The number of models for the minority and majority classifier are not the same.
         Please check your models within the minorityDir and majorityDir that the
         same seeds have been used for the generation of a minority and a majority classifier.")
  } else if (!all.equal(selectedDirsMajority, selectedDirsMinority) ) {
    base::stop("Please make sure you run the crossvalidation with the same seed for complementary classifications,
         and store them in the same directory.")
  }
  } else {
    selectedDirsMajority <- majorityDir
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
        base::print("Either change the column with the tumor subtype labels to the name: ", classColumn)
        base::stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
      } else if (higherClassColumn %notin% base::colnames(metaDataTest)) {

        base::print("Please note that the wanted column for the tumor type labels cannot be found within 'metaDataTest'.")
        base::print("Either change the column with the tumor type labels to the name: ", higherClassColumn)
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
    if (crossValidation == F & subtype == F) {
      predictionsMMFinal$originalCall <- metaDataTest[base::rownames(predictionsMMFinal), higherClassColumn]
      } else if (crossValidation == F) {
      predictionsMMFinal$originalCall <- metaDataTest[base::rownames(predictionsMMFinal), classColumn]
      }


    if (subtype == T) {
      fractionsCorrect <- getAccuraciesPerTumorTypeSize(predictionsMMFinal,
                                                        metaDataRef = minority$metaDataRef,
                                                        classColumn = classColumn,
                                                        rounding = rounding,
                                                        probabilityThreshold = probabilityThreshold)
    } else {
      fractionsCorrect <- getAccuraciesPerTumorTypeSize(predictionsMMFinal,
                                                        metaDataRef = minority$metaDataRef,
                                                        classColumn = higherClassColumn,
                                                        rounding = rounding,
                                                        probabilityThreshold = probabilityThreshold)

    }
    fractionsCorrect$seed <- i
    if (i == 1) {
      accuracyDF <- fractionsCorrect
    } else {
      accuracyDF <- base::rbind(accuracyDF, fractionsCorrect)
    }
  }

  meanNumbers <- accuracyDF %>%
    dplyr::group_by(nCases) %>%
    dplyr::summarise(
      meanFractionCorrect = base::mean(fractionCorrect),
      meanFractionCorrectFiltered = base::mean(fractionCorrectFiltered),
      meanFractionIncorrect = base::mean(1 - fractionCorrect),
      meanFractionIncorrectFiltered = base::mean(1 - fractionCorrectFiltered),
      sdFractionCorrect = stats::sd(fractionCorrect),
      sdFractionCorrectFiltered = stats::sd(fractionCorrectFiltered),
      meanCasesFiltered = base::round(base::mean(nSamplesFiltered), digits = 0),
      meanPrecision = base::mean(Precision),
      meanF1 = base::mean(F1),
      meanFractionCorrect2 = base::mean(fractionCorrect2),
      meanFractionCorrect3 = base::mean(fractionCorrect3),
      sdFractionCorrect2 = stats::sd(fractionCorrect2),
      sdFractionCorrect3 = stats::sd(fractionCorrect3),
      meanRecall = base::mean(Recall),
      medianF1 = stats::median(F1),
      sdPrecision = stats::sd(Precision),
      sdF1 = stats::sd(F1),
      sdRecall = stats::sd(Recall),
      meanSamples = base::mean(nSamples)
    )

  if (crossValidation == T) {
    meanNumbers$type <- "Train"
  } else {
    meanNumbers$type <- "Test"
  }

  meanNumbers$subtype <- subtype

  return(meanNumbers)
}
