#' Calculate the influence of tumor heterogeneity on performance M&M
#'
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param trainOrTest Do you want to calculate for the reference cohort ("Train") or independent test set ("Test")?
#' @param nSeeds How many seeds was the cross-validation setup run with?
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param metaDataTest Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the test set.
#' @param probabilityScoreTumor What is the probability score threshold you would like to use to call a classification 'confident' on the tumor type level?
#' @param probabilityScoreSubtype What is the probability score threshold you would like to use to call a classification 'confident' on the tumor subtype level?
#' @param metaDataFFPE Metadata file containing the links between the patients and the tumor (sub)type diagnosis for the FFPE samples.
#' @param predictionsFFPE Predictions for the FFPE samples on the tumor type level.
#' @param predictionsFFPESubtype Predictions for the FFPE samples on the tumor subtype level.
#'
#' @return Dataframe containing the results concerning the samples from different sources of tumor heterogeneity for all train or test datasets.
#' @export
#' @import magrittr
#'
calculateMeanAndSDInfluenceTumorHeterogeneity <- function(
    minorityDir,
    majorityDir,
    metaDataRef,
    metaDataTest = NA,
    probabilityScoreTumor,
    probabilityScoreSubtype,
    metaDataFFPE = NA,
    predictionsFFPE,
    predictionsFFPESubtype
  ) {


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
    trainOrTest <- "Train"
  } else {
    crossValidation <- F
    trainOrTest <- "Test"
  }

  for (i in seq(1:length(selectedDirsMajority))) {
    if (trainOrTest == "Train") {
      #crossValidation <- T
      minorityDoc <- paste0(minorityDir, "seed",i, "/crossValidationMinorityResults.rds")
      majorityDoc <- paste0(majorityDir, "seed",i, "/crossValidationMajorityResults.rds")
    } else {
      #crossValidation <- F
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
          print("Either change the column with the tumor subtype labels to the name: ", classColumn)
          stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
        } else if (higherClassColumn %notin% colnames(metaDataTest)) {

          print("Please note that the wanted column for the tumor type labels cannot be found within 'metaDataTest'.")
          print("Either change the column with the tumor type labels to the name: ", higherClassColumn)
          stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
        } else {
          print(paste0("Found columns ", classColumn, " and ", higherClassColumn, " within metaDataTest specifying the tumor subtype, and tumor type."))

          print("No original call found, adding it from metaDataTest")
        }
      }
    }

    if (trainOrTest == "Train") {
      crossValidation <- T
    } else {
      crossValidation <- F
    }
    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          subtype = F
    )

    predictionsMM <- predictionsMMFinalList$predictionsMMFinal

    #predictionsMM %<>% filter(rownames(.) %notin% throwOut)


    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          subtype = T
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

     statusDFLonger <- InfluenceTumorHeterogeneity(metaDataRef = metaDataTest,
                                                   predictionsMM = predictionsMM,
                                                   predictionsMMSubtype = predictionsMMSubtype,
                                                   trainOrTest = trainOrTest,
                                                   probabilityScoreTumor = probabilityScoreTumor,
                                                   probabilityScoreSubtype = probabilityScoreSubtype,
                                                   metaDataFFPE = metaDataFFPE,
                                                   predictionsFFPE = predictionsFFPE,
                                                   predictionsFFPESubtype = predictionsFFPESubtype)
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
      sdAccuracy = sd(accuracy),
      sdPrecision = sd(precision),
      sdRecall = sd(recall)
    )

  return(meanNumbers)

}
