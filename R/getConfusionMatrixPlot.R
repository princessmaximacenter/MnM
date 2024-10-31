#' Get confusion matrix information
#'
#' This function is designed to extract which tiles of the confusion matrix
#' do contain reference-prediction information, and which do not.

#' IMPORTANT: The order of the tumor types within the plot will be equal to the order of the
#' tumors within the _abbreviations_ dataframe. Make sure you specify the
#' order of the tumor types well within the abbreviations, not only concerning the domains,
#' but also the subsequent tumor types order and tumor subtypes order.
#'
#' @param minorityDir Directory in which the Minority classifier model(s) are stored.
#' @param majorityDir in which the Majority classifier model(s) are stored.
#' @param abbreviations Dataframe containing the links between the tumor subtype and their abbreviations ($abbreviationSubtype),
#' the tumor types and their abbreviations ($abbreviationTumorType), and the domain.
#' @param probabilityThreshold Which probability score do you want to use as a cutoff for calling classifications confident or not?
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' @import caret
#' @return List containing all tiles that contain classifications with their frequencies ($confusionPlotDF),
#' which tiles do not contain classifications ($nonAvailableTiles),
#' what were the abbreviations used as input for generating the  ($abbreviations),
#' what is the column in the abbreviations that contains the tumor subtype labels ($classColumn),
#' what is the column in the abbreviations that contains the tumor type labels ($higherClassColumn),
#' and lastly whether the tiles with classification frequencies were based on tumor type($subtype = F) or subtype level ($subtype = T).

#' @export
#'
getConfusionMatrixPlot <- function(minorityDir,
                                   majorityDir,
         abbreviations = NA,
         subtype,
         probabilityThreshold) {

  `%notin%` <- base::Negate(`%in%`)


  predictionsMMAverageList <- combineSeedPredictions(
    minorityDir = minorityDir,
    majorityDir = majorityDir,
    subtype = subtype
  )

  metaDataRef <- predictionsMMAverageList$metaDataRef
  classColumn <- predictionsMMAverageList$metaDataRun$classColumn
  higherClassColumn <- predictionsMMAverageList$metaDataRun$higherClassColumn
  domainColumn <- predictionsMMAverageList$metaDataRun$domainColumn

  predictionsMM <- predictionsMMAverageList$predictionsMMFinal

  if (base::is.na(abbreviations)[1] & !base::is.na(domainColumn)[1]) {

    abbreviations <- metaDataRef[, c(domainColumn, classColumn, higherClassColumn)] %>% base::unique()
    abbreviations$abbreviationSubtype <- abbreviations[,classColumn]
    abbreviations$abbreviationTumorType <- abbreviations[,higherClassColumn]

    base::print("You have not supplied any abbreviations (abbreviations). If you would like to use abbreviations,
          please generate a dataframe with the following columns and abbreviations within abbreviationSubtype and abbreviationTumorType:")
    base::print(abbreviations[1:4,])
  } else if (base::is.na(abbreviations)[1] ) {
    abbreviations <- metaDataRef[, c(classColumn, higherClassColumn)] %>% base::unique()
    abbreviations$domainColumn <- " "
    abbreviations$abbreviationSubtype <- abbreviations[,classColumn]
    abbreviations$abbreviationTumorType <- abbreviations[,higherClassColumn]

    base::print("You have not supplied any abbreviations (abbreviations). If you would like to use abbreviations,
          please generate a dataframe with the following columns and abbreviations within abbreviationSubtype and abbreviationTumorType:")
    base::print(abbreviations[1:4,])
  }


  if (base::is.na(domainColumn)[1]) {
    domainColumn <- "domainColumn"

    metaDataRef$domainColumn <- " "
  }

  if (subtype == T) {
    abbreviations$abbreviation <- abbreviations$abbreviationSubtype
    chosenClassColumn <- classColumn


  } else {
    abbreviations$abbreviation <- abbreviations$abbreviationTumorType
    chosenClassColumn <- higherClassColumn
  }

  if ("Low confidence" %notin% abbreviations[,"abbreviation"]) {
    notClassifiedLine <- data.frame(classColumn = "Low confidence",
                                    abbreviationSubtype = NA,
                                    abbreviationTumorType = NA,
                                    higherClassColumn = "Low confidence",
                                    domainColumn = NA,
                                    abbreviation = "Low confidence")

    colnames(notClassifiedLine) <- base::gsub("classColumn", classColumn, base::colnames(notClassifiedLine))
    colnames(notClassifiedLine) <- base::gsub("higherClassColumn", higherClassColumn, base::colnames(notClassifiedLine))
    colnames(notClassifiedLine) <- base::gsub("domainColumn", domainColumn, base::colnames(notClassifiedLine))

    abbreviations <- base::rbind(notClassifiedLine, abbreviations[,base::colnames(notClassifiedLine)])

  }

  predictionsMMFiltered <- predictionsMM %>% dplyr::filter(probability1 > probabilityThreshold)

  tumorConfusionMatrix <- caret::confusionMatrix(factor(predictionsMMFiltered$predict,
                                                 levels = base::unique(c(predictionsMM$originalCall, predictionsMM$predict))),
                                          factor(predictionsMMFiltered$originalCall,
                                                 levels = base::unique(c(predictionsMM$originalCall,
                                                                                                       predictionsMM$predict))),
                                          dnn = c("Prediction", "Reference"))
  predictionFrequencies <- tumorConfusionMatrix$table %>% base::as.data.frame()
  predictionFrequencies$Prediction <- base::as.character(predictionFrequencies$Prediction)
  predictionFrequencies$Reference <- base::as.character(predictionFrequencies$Reference)
  missingTumors <- predictionFrequencies %>% dplyr::filter(Reference == Prediction,
                                                    Freq == 0)

  predictionFrequencies %<>% dplyr::filter(Freq != 0)

  for (j in base::seq(1:base::nrow(missingTumors))) {


    if (missingTumors$Reference[j] %notin% base::unique(predictionsMMFiltered$originalCall)) {
      missingTumors <- base::rbind(missingTumors, base::data.frame(Prediction = "Low confidence",
                                                       Reference = missingTumors$Reference[j],
                                                       Freq = (base::nrow(predictionsMM[predictionsMM$originalCall == missingTumors$Reference[j], ]) -
                                                                 base::nrow(predictionsMMFiltered[predictionsMMFiltered$originalCall == missingTumors$Reference[j], ])))
      )

    }
  }

  for (j in base::seq(1:base::nrow(missingTumors))) {
    missingTumors$Prediction[j] <- abbreviations[abbreviations[,chosenClassColumn] == missingTumors$Prediction[j], "abbreviation", drop = T] %>% base::unique()
    missingTumors$Reference[j] <- abbreviations[abbreviations[,chosenClassColumn] == missingTumors$Reference[j], "abbreviation", drop = T] %>% base::unique()
  }




  difPredictions <- base::unique(predictionFrequencies$Reference) %>% base::as.character()


  linkClassAndDomain <- metaDataRef[ , c( chosenClassColumn, domainColumn)] %>% base::unique()

  for (i in base::seq(1:base::length(difPredictions))) {
    total <- predictionFrequencies %>% dplyr::filter(Reference == difPredictions[i])
    total$Domain <- linkClassAndDomain[linkClassAndDomain[,chosenClassColumn] == total$Reference[1], domainColumn]
    totalNum <- base::sum(total$Freq)
    metaTotal <-predictionsMM %>%
      dplyr::filter(originalCall == difPredictions[i]) %>%
      base::nrow()
    notClassified <- metaTotal - totalNum


    for (j in base::seq(1:base::nrow(total))) {
      total$Prediction[j] <- abbreviations[abbreviations[,chosenClassColumn] == total$Prediction[j], "abbreviation", drop = T] %>% base::unique()
      total$Reference[j] <- abbreviations[abbreviations[,chosenClassColumn] == total$Reference[j], "abbreviation", drop = T] %>% base::unique()
    }

    newLine <- data.frame(Prediction = "Low confidence",
                          Reference = base::as.character(total$Reference[1]),
                          Freq = notClassified,
                          Domain = base::unique(total$Domain)
    )

    total <- base::rbind(total, newLine)

    if (i == 1) {
      confusionPlotDF <- total
    } else {
      confusionPlotDF <- base::rbind(confusionPlotDF, total)
    }
  }

  linkClassAndDomain$abbreviation <- NA

  for (i in base::seq(1:base::nrow(linkClassAndDomain))) {
    if (linkClassAndDomain[i, chosenClassColumn] %in% abbreviations[,chosenClassColumn]) {
      linkClassAndDomain[i, "abbreviation"] <-
        abbreviations[abbreviations[,chosenClassColumn] == linkClassAndDomain[i,chosenClassColumn], "abbreviation"] %>%
        base::unique()
    }
  }


  confusionPlotDF$Prediction <- base::factor(confusionPlotDF$Prediction,
                                             levels = base::unique(abbreviations$abbreviation[abbreviations$abbreviation %in%
                                                                                                base::unique(c(confusionPlotDF$Reference, confusionPlotDF$Prediction,
                                                                                                    missingTumors$Reference))]))

  missingNotClassified <- base::levels(confusionPlotDF$Prediction)[base::levels(confusionPlotDF$Prediction) %notin% confusionPlotDF[confusionPlotDF$Prediction == "Low confidence","Reference"] ]
  missingNotClassified <- missingNotClassified[missingNotClassified != "Low confidence"]

  for (i in base::seq(1:base::length(missingNotClassified))) {
    missingNotClassifiedDF <- base::data.frame(Prediction = "Low confidence",
                                         Reference = missingNotClassified[i],
                                         Freq = 0,
                                         Domain = linkClassAndDomain[linkClassAndDomain[,"abbreviation"] == missingNotClassified[i],domainColumn]    )
    if (i == 1)  {
      missingClassifiedDF <- missingNotClassifiedDF
    } else {
      missingClassifiedDF <- base::rbind(missingClassifiedDF,
                                   missingNotClassifiedDF)
    }
  }

  if (subtype == T) {
    confusionPlotDF$Reference <- base::factor(confusionPlotDF$Reference,
                                 levels = base::levels(confusionPlotDF$Prediction))
  } else {
    confusionPlotDF$Reference <- base::factor(confusionPlotDF$Reference,
                                 levels = base::levels(confusionPlotDF$Prediction)[-1])
  }

  missingTumors$Domain <- NA
  for (i in base::seq(1:base::nrow(missingTumors))) {
    missingTumors$Domain[i] = linkClassAndDomain[linkClassAndDomain$abbreviation == missingTumors$Reference[i], domainColumn]
  }
  missingTumors$Prediction <- base::factor(missingTumors$Prediction, levels = base::levels(confusionPlotDF$Prediction))
  missingTumors$Reference <- base::factor(missingTumors$Reference, levels = base::levels(confusionPlotDF$Reference))

  confusionPlotDF <- base::rbind(confusionPlotDF, missingTumors)


    nonAvailableTiles <- getNonAvailableTiles(predictionsMM = predictionsMM,
                                              abbreviations = abbreviations,
                                              classColumn = chosenClassColumn,
                                              probabilityThreshold = probabilityThreshold)

  confusionPlotInfo <- base::list(confusionPlotDF = confusionPlotDF,
                            nonAvailableTiles = nonAvailableTiles,
                            abbreviations = abbreviations,
                            classColumn = classColumn,
                            higherClassColumn = higherClassColumn,
                            subtype = subtype
                            )


  return(confusionPlotInfo)
}
