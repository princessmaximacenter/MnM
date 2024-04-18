#' Calculate performance of M&M on total train and test set
#'
#' @param classColumn Column in the original metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the original metadata file that contains the tumor type labels.
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param metaDataTest Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the test set.
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' @param rounding Do you want rounded numbers for the performance scores?
#' @param probabilityThreshold What is the probability score threshold you would like to use to call a classification 'confident'?
#'
#' @return Dataframe containing the average performance of tumor classifications within a certain frequency range (nCases).
#' Included are the averages for the percentage of correctly and incorrectly classified samples ($meanFractionCorrect and $meanFractionIncorrect),
#' correctly and incorrectly classified 'confident' samples ($meanFractionCorrectFiltered and $meanFractionIncorrectFiltered),
#' averaged precision values for the average of all tumor entities within the frequency range ($meanPrecision),
#' averaged recall values for the average of all tumor entities within the frequency range ($meanRecall),
#' and averaged F1 scores for the average of all tumor entities within the frequency range ($meanF1).
#' Please note that the precision, F1 and recall are calculated for the confident sample classifications only.
#'The total amount of samples within each frequency range ($meanSamples) is also specified.
#' @export
#' @import magrittr dplyr
calculateMeanAndSDAccuracy <- function(classColumn,
                                       higherClassColumn,
                                       minorityDir,
                                       majorityDir,
                                       metaDataTest = NA,
         subtype = F,
         rounding = F,
         probabilityThreshold
         ) {

  `%notin%` <- Negate(`%in%`)
  if ("minorityClassifierResult.rds" %notin% list.files(minorityDir)) {
    crossValidation <- T
  allDirsMinority <- list.dirs(minorityDir, recursive = F)
  allDirsMajority <- list.dirs(majorityDir, recursive = F)
  selectedDirsMinority <- allDirsMinority[grep("seed", allDirsMinority)]
  selectedDirsMajority <- allDirsMajority[grep("seed", allDirsMajority)]

  print(paste0("Found ",length(selectedDirsMajority), " directories with different cross-validation runs.",
  " Calculating average performance values for all combined."))
  if (length(selectedDirsMinority) != length(selectedDirsMajority)) {
    stop("The number of models for the minority and majority classifier are not the same.
         Please check your models within the minorityDir and majorityDir that the
         same seeds have been used for the generation of a minority and a majority classifier.")
  } else if (!all.equal(selectedDirsMajority, selectedDirsMinority) ) {
    stop("Please make sure you run the crossvalidation with the same seed for complementary classifications,
         and store them in the same directory.")
  }
  } else {
    selectedDirsMajority <- majorityDir
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

    predictionsMMFinalList <- integrateMM(minority = minority,
                                      majority = majority,
                                      subtype = subtype,
                                      classColumn = classColumn,
                                      higherClassColumn = higherClassColumn)

    predictionsMMFinal <- predictionsMMFinalList$predictionsMMFinal
    if (crossValidation == F & subtype == F) {
      predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), higherClassColumn]
      } else if (crossValidation == F) {
      predictionsMMFinal$originalCall <- metaDataTest[rownames(predictionsMMFinal), classColumn]
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
      accuracyDF <- rbind(accuracyDF, fractionsCorrect)
    }
  }

  meanNumbers <- accuracyDF %>% group_by(nCases) %>%
    summarise(
      meanFractionCorrect = mean(fractionCorrect),
      meanFractionCorrectFiltered = mean(fractionCorrectFiltered),
      meanFractionIncorrect = mean(1 - fractionCorrect),
      meanFractionIncorrectFiltered = mean(1 - fractionCorrectFiltered),
      sdFractionCorrect = sd(fractionCorrect),
      sdFractionCorrectFiltered = sd(fractionCorrectFiltered),
      meanCasesFiltered = round(mean(nSamplesFiltered), digits = 0),
      meanPrecision = mean(Precision),
      meanF1 = mean(F1),
      meanFractionCorrect2 = mean(fractionCorrect2),
      meanFractionCorrect3 = mean(fractionCorrect3),
      sdFractionCorrect2 = sd(fractionCorrect2),
      sdFractionCorrect3 = sd(fractionCorrect3),
      meanRecall = mean(Recall),
      medianF1 = median(F1),
      sdPrecision = sd(Precision),
      sdF1 = sd(F1),
      sdRecall = sd(Recall),
      meanSamples = mean(nSamples)
    )

  return(meanNumbers)
}
