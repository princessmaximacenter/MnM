#' Title
#'
#' @param metaDataRef
#' @param predictionsMM
#' @param predictionsMMSubtype
#' @param trainOrTest
#' @param probabilityScoreTumor
#' @param probabilityScoreSubtype
#' @param metaDataFFPE
#' @param predictionsFFPE
#' @param predictionsFFPESubtype
#'
#' @return
#' @export
#'
#' @examples
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
