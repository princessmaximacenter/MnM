InfluenceTumorHeterogeneity <- function(metaDataRef,
                                        predictionsMM,
                                        predictionsMMSubtype,
                                        trainOrTest,
                                       # predictionsMMTest,
                                       # predictionsMMTestSubtype,
                                        probabilityScoreTumor,
                                        probabilityScoreSubtype,
                                       testSet) {

  tumorStatusTrain <- primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataRef,
                                                         predictionsMM = predictionsMM,
                                                         columnOfInterest = "Status",
                                                         labelTypes = c("Primary", "Recurrence", "Metastasis"),
                                                         probabilityScore = probabilityScoreTumor,
                                                         includeAllRow = T)
#  tumorStatusTrain$trainOrTest <- trainOrTest
  tumorStatusTrain$type <- "TumorType"

  treatmentStatusTrain <- primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataRef,
                                                             predictionsMM = predictionsMM,
                                                             columnOfInterest = "SystemicTreatment",
                                                             labelTypes = c("Yes", "No"),
                                                             probabilityScore = probabilityScoreTumor)
  # if (crossValidation == T) {
  # treatmentStatusTrain$trainOrTest <- "Train"
  # } else {
  #   treatmentStatusTrain$trainOrTest <- "Test"
  # }

  treatmentStatusTrain$type <- "TumorType"
  treatmentStatusTrain$labelType <- c("Systemic treatment",
                                      "No treatment")


  # metaDataTest <- newTestSet$metaData %>% filter(rownames(.) %notin% throwOut)
  # tumorStatusTest <- primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataTest,
  #                                                       predictionsMM = predictionsMMTestFiltered,
  #                                                       probabilityScore = probabilityScoreTumor,
  #                                                       columnOfInterest = "Status",
  #                                                       labelTypes = c("Primary", "Recurrence", "Metastasis")
  # )


 # tumorStatusTest$trainOrTest <- "Test"
 # tumorStatusTest$type <- "TumorType"



  # treatmentStatusTest <- primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataTest,
  #                                                           predictionsMM = predictionsMMTestFiltered,
  #                                                           probabilityScore = probabilityScoreTumor,
  #                                                           columnOfInterest = "SystemicTreatment",
  #                                                           labelTypes = c("Yes", "No")
  # )
  # treatmentStatusTest$trainOrTest <- "Test"
  # treatmentStatusTest$type <- "TumorType"
  # treatmentStatusTest$labelType <- c("Systemic treatment",
  #                                    "No treatment")

  tumorStatusSubtypeTrain <-  primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataRef,
                                                                 predictionsMM = predictionsMMSubtype,
                                                                 probabilityScore = probabilityScoreSubtype,
                                                                 columnOfInterest = "Status",
                                                                 labelTypes = c("Primary", "Recurrence", "Metastasis"),
                                                                 includeAllRow = T
  )
  tumorStatusSubtypeTrain$type <- "TumorSubtype"
  #tumorStatusSubtypeTrain$trainOrTest <- trainOrTest


  treatmentStatusSubtypeTrain <-  primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataRef,
                                                                     predictionsMM = predictionsMMSubtype,
                                                                     probabilityScore = probabilityScoreSubtype,
                                                                     columnOfInterest = "SystemicTreatment",
                                                                     labelTypes = c("Yes", "No")
  )
  treatmentStatusSubtypeTrain$type <- "TumorSubtype"
  #treatmentStatusSubtypeTrain$trainOrTest <- "Train"
  treatmentStatusSubtypeTrain$labelType <- c("Systemic treatment",
                                             "No treatment")
  if (trainOrTest == "Test") {
  predictionsMMTestPreviously <- read.table("~/Documents/data/tumorClassification/output/finalResults/MnMResults/predictionsMMTest.tsv",
                                            sep = "\t",
                                            header = T)

  predictionsMMSubtypeTestPreviously <- read.table("~/Documents/data/tumorClassification/output/finalResults/MnMResults/predictionsMMSubtypeTest.tsv",
                                                   sep = "\t",
                                                   header = T)

  FFPEResults <- primaryRecurrenceMetastasisNumbers(metaDataRef = testSet$metaData,
                                                    predictionsMM = predictionsMMTestPreviously,
                                                    labelTypes = "FFPE",
                                                    columnOfInterest = "Disease_sub_specification2",
                                                    probabilityScore = 0.8)



  FFPEResultsSubtype <- primaryRecurrenceMetastasisNumbers(metaDataRef = testSet$metaData,
                                                           predictionsMM = predictionsMMSubtypeTestPreviously,
                                                           labelTypes = "FFPE",
                                                           columnOfInterest = "Disease_sub_specification2",
                                                           probabilityScore = 0.7)

  FFPEResults$type <- "TumorType"





  FFPEResultsSubtype$type <- "TumorSubtype"

  }
  # tumorStatusSubtypeTest <-  primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataTest,
  #                                                               predictionsMM = predictionsMMSubtypeTestFiltered,
  #                                                               probabilityScore = probabilityScoreSubtype,
  #                                                               columnOfInterest = "Status",
  #                                                               labelTypes = c("Primary", "Recurrence", "Metastasis")
  # )
  # tumorStatusSubtypeTest$type <- "TumorSubtype"
  # tumorStatusSubtypeTest$trainOrTest <- "Test"


  # treatmentStatusSubtypeTest <-  primaryRecurrenceMetastasisNumbers(metaDataRef = metaDataTest,
  #                                                                   predictionsMM = predictionsMMSubtypeTestFiltered,
  #                                                                   probabilityScore = probabilityScoreSubtype,
  #                                                                   columnOfInterest = "SystemicTreatment",
  #                                                                   labelTypes = c("Yes", "No")
  # )
 # treatmentStatusSubtypeTest$type <- "TumorSubtype"
 # treatmentStatusSubtypeTest$trainOrTest <- "Test"
 # treatmentStatusSubtypeTest$labelType <- c("Systemic treatment",
 #                                           "No treatment")


  statusDF <- rbind(tumorStatusTrain,
                   # tumorStatusTest,
                    tumorStatusSubtypeTrain,
                   # tumorStatusSubtypeTest,
                    treatmentStatusTrain,
                    #treatmentStatusTest,
                    treatmentStatusSubtypeTrain
                   # treatmentStatusSubtypeTest
  )

  if (trainOrTest == "Test") {
    statusDF <- rbind(statusDF, FFPEResults, FFPEResultsSubtype)

  }
  statusDF %<>% unique()

  #statusDFLonger <- statusDF %>% pivot_longer(cols = c("accuracy",
    #                                                   "precision",
     #                                                  "recall"
 # ),
  #names_to = "measurementType"
 # )
 # statusDFLonger$valuePercent <- paste0(round(statusDFLonger$value * 100,2), "%")
  # statusDFLonger$labelType <- factor(statusDFLonger$labelType,
  #                                    levels = c("Primary",
  #                                               "Recurrence",
  #                                               "Metastasis",
  #                                               "No treatment",
  #                                               "Systemic treatment"))

  statusDF$labelType <- factor(statusDF$labelType,
                                     levels = c("All", "Primary",
                                                "Recurrence",
                                                "Metastasis",
                                                "No treatment",
                                                "Systemic treatment",
                                                "FFPE"))

  #statusDFLonger$trainOrTest <- trainOrTest
  statusDF$trainOrTest <- trainOrTest
 # statusDFLonger$trainOrTest <- factor(statusDFLonger$trainOrTest,
  #                                     levels = c("Train", "Test"))

  # statusDFLonger$type <- factor(statusDFLonger$type,
  #                               levels = c("TumorType", "TumorSubtype"))

  statusDF$type <- factor(statusDF$type,
                          levels = c("TumorType", "TumorSubtype"))
  return(statusDF)

}
