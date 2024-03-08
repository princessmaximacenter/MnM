#' Calculate the influence of tumor heterogeneity on performance M&M
#'
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param trainOrTest Do you want to calculate for the reference cohort ("Train") or independent test set ("Test")?
#' @param nSeeds How many seeds was the cross-validation setup run with?
#' @param nModels How many models were created for the majority voting system?
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param metaDataTest Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the test set.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param probabilityScoreTumor What is the probability score threshold you would like to use to call a classification 'confident' on the tumor type level?
#' @param probabilityScoreSubtype What is the probability score threshold you would like to use to call a classification 'confident' on the tumor subtype level?
#' @param metaDataFFPE Metadata file containing the links between the patients and the tumor (sub)type diagnosis for the FFPE samples.
#' @param predictionsFFPE Predictions for the FFPE samples on the tumor type level.
#' @param predictionsFFPESubtype Predictions for the FFPE samples on the tumor subtype level.
#' @param throwOut Are there samples you would like to remove from the test set due to poor data quality? If so, add their rownames here.
#'
#' @return Dataframe containing the results concerning the samples from different sources of tumor heterogeneity for all train or test datasets.
#' @export
#' @import magrittr
#'
calculateMeanAndSDInfluenceTumorHeterogeneity <- function(
    minorityDir,
    majorityDir,
    trainOrTest,
    nSeeds,
    nModels,
    metaDataRef,
    metaDataTest = NA,
    classColumn,
    higherClassColumn,
    probabilityScoreTumor,
    probabilityScoreSubtype,
    throwOut,
    metaDataFFPE = NA,
    predictionsFFPE,
    predictionsFFPESubtype
  ) {


  for (i in seq(1:nSeeds)) {
    if (trainOrTest == "Train") {
      crossValidation <- T
      minorityDoc <- paste0(minorityDir, "seed",i, "/crossValidationMinorityResults.rds")
      majorityDoc <- paste0(majorityDir, "seed",i, "/crossValidationMajorityResults.rds")
    } else {
      crossValidation <- F
      minorityDoc <- minorityDir
      majorityDoc <- majorityDir
    }


    minority <- readRDS(minorityDoc)
    majority <- readRDS(majorityDoc)

    if (trainOrTest == "Train") {
      crossValidation <- T
    } else {
      crossValidation <- F
    }
    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          nModels = nModels,
                                          subtype = F,
                                          metaDataRef = metaDataRef,
                                          classColumn = classColumn,
                                          higherClassColumn = higherClassColumn,
                                          crossValidation = crossValidation
    )

    predictionsMM <- predictionsMMFinalList$predictionsMMFinal

    #predictionsMM %<>% filter(rownames(.) %notin% throwOut)


    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          nModels = nModels,
                                          subtype = T,
                                          metaDataRef = metaDataRef,
                                          classColumn = classColumn,
                                          higherClassColumn = higherClassColumn,
                                          crossValidation = crossValidation
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
     predictionsMMSubtype %<>% filter(rownames(.) %notin% throwOut)
     predictionsMM %<>% filter(rownames(.) %notin% throwOut)
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
