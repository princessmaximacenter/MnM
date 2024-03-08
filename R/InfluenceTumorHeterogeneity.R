#'  Calculate results for single dataset concerning tumor heterogeneity
#'
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the reference cohort.
#' @param predictionsMM Dataframe containing the top 3 final tumor type classification labels ($predict{2,3}) with their accompanying probability scores
#' ($probability{1,2,3} and the original diagnosis label ($originalCall).
#' @param predictionsMMSubtype Dataframe containing the top 3 final tumor subtype classification labels ($predict{2,3}) with their accompanying probability scores
#' ($probability{1,2,3} and the original diagnosis label ($originalCall).
#' @param trainOrTest Do you want to calculate for the reference cohort ("Train") or independent test set ("Test")?
#' @param probabilityScoreTumor What is the probability score threshold you would like to use to call a classification 'confident' on the tumor type level?
#' @param probabilityScoreSubtype What is the probability score threshold you would like to use to call a classification 'confident' on the tumor subtype level?
#' @param metaDataFFPE Metadata file containing the links between the patients and the tumor (sub)type diagnosis for the FFPE samples.
#' @param predictionsFFPE Predictions for the FFPE samples on the tumor type level.
#' @param predictionsFFPESubtype Predictions for the FFPE samples on the tumor subtype level.
#'
#' @return Dataframe containing the results for different sources of heterogeneity within one dataset.
#'
InfluenceTumorHeterogeneity <- function(metaDataRef,
                                        predictionsMM,
                                        predictionsMMSubtype,
                                        trainOrTest,
                                        probabilityScoreTumor,
                                        probabilityScoreSubtype,
                                       metaDataFFPE = NA,
                                       predictionsFFPE,
                                       predictionsFFPESubtype
                                       ) {

  tumorStatusTrain <- primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataRef,
                                                         predictionsMM = predictionsMM,
                                                         columnOfInterest = "Status",
                                                         labelTypes = c("Primary", "Recurrence", "Metastasis"),
                                                         probabilityScore = probabilityScoreTumor,
                                                         includeAllRow = T)

  tumorStatusTrain$type <- "TumorType"

  treatmentStatusTrain <- primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataRef,
                                                             predictionsMM = predictionsMM,
                                                             columnOfInterest = "SystemicTreatment",
                                                             labelTypes = c("Yes", "No"),
                                                             probabilityScore = probabilityScoreTumor)


  treatmentStatusTrain$type <- "TumorType"
  treatmentStatusTrain$labelType <- c("Systemic treatment",
                                      "No treatment")

  tumorStatusSubtypeTrain <-  primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataRef,
                                                                 predictionsMM = predictionsMMSubtype,
                                                                 probabilityScore = probabilityScoreSubtype,
                                                                 columnOfInterest = "Status",
                                                                 labelTypes = c("Primary", "Recurrence", "Metastasis"),
                                                                 includeAllRow = T
  )
  tumorStatusSubtypeTrain$type <- "TumorSubtype"


  treatmentStatusSubtypeTrain <-  primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataRef,
                                                                     predictionsMM = predictionsMMSubtype,
                                                                     probabilityScore = probabilityScoreSubtype,
                                                                     columnOfInterest = "SystemicTreatment",
                                                                     labelTypes = c("Yes", "No")
  )
  treatmentStatusSubtypeTrain$type <- "TumorSubtype"
  treatmentStatusSubtypeTrain$labelType <- c("Systemic treatment",
                                             "No treatment")

  statusDF <- rbind(tumorStatusTrain,
                    tumorStatusSubtypeTrain,
                    treatmentStatusTrain,
                    treatmentStatusSubtypeTrain)


  if (!is.na(metaDataFFPE)[1]) {

  FFPEResults <- primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataFFPE,
                                                    predictionsMM = predictionsFFPE,
                                                    labelTypes = "FFPE",
                                                    columnOfInterest = "general_observation",
                                                    probabilityScore = probabilityScoreTumor)



  FFPEResultsSubtype <- primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataFFPE,
                                                           predictionsMM = predictionsFFPESubtype,
                                                           labelTypes = "FFPE",
                                                           columnOfInterest = "general_observation",
                                                           probabilityScore = probabilityScoreSubtype)

  FFPEResults$type <- "TumorType"





  FFPEResultsSubtype$type <- "TumorSubtype"


  statusDF <- rbind(statusDF,
                    FFPEResults,
                    FFPEResultsSubtype)
  }
  statusDF %<>% unique()


  statusDF$labelType <- factor(statusDF$labelType,
                                     levels = c("All", "Primary",
                                                "Recurrence",
                                                "Metastasis",
                                                "No treatment",
                                                "Systemic treatment",
                                                "FFPE"))

  statusDF$trainOrTest <- trainOrTest
  statusDF$type <- factor(statusDF$type,
                          levels = c("TumorType", "TumorSubtype"))
  return(statusDF)

}
