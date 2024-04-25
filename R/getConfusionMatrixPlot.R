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
#' @param majorityDirDirectory in which the Majority classifier model(s) are stored.
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

  `%notin%` <- Negate(`%in%`)


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

  if (is.na(abbreviations)[1] & !is.na(domainColumn)[1]) {

    abbreviations <- metaDataRef[, c(domainColumn, classColumn, higherClassColumn)] %>% base::unique()
    abbreviations$abbreviationSubtype <- abbreviations[,classColumn]
    abbreviations$abbreviationTumorType <- abbreviations[,higherClassColumn]

    print("You have not supplied any abbreviations (abbreviations). If you would like to use abbreviations,
          please generate a dataframe with the following columns and abbreviations within abbreviationSubtype and abbreviationTumorType:")
    print(abbreviations[1:4,])
  } else if (is.na(abbreviations)[1] ) {
    abbreviations <- metaDataRef[, c(classColumn, higherClassColumn)] %>% base::unique()
    abbreviations$domainColumn <- " "
    abbreviations$abbreviationSubtype <- abbreviations[,classColumn]
    abbreviations$abbreviationTumorType <- abbreviations[,higherClassColumn]

    print("You have not supplied any abbreviations (abbreviations). If you would like to use abbreviations,
          please generate a dataframe with the following columns and abbreviations within abbreviationSubtype and abbreviationTumorType:")
    print(abbreviations[1:4,])
  }


  if (is.na(domainColumn)[1]) {
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

  if ("Not classified" %notin% abbreviations[,"abbreviation"]) {
    notClassifiedLine <- data.frame(classColumn = "Not classified",
                                    abbreviationSubtype = NA,
                                    abbreviationTumorType = NA,
                                    higherClassColumn = "Not classified",
                                    domainColumn = NA,
                                    abbreviation = "Not classified")

    colnames(notClassifiedLine) <- gsub("classColumn", classColumn, colnames(notClassifiedLine))
    colnames(notClassifiedLine) <- gsub("higherClassColumn", higherClassColumn, colnames(notClassifiedLine))
    colnames(notClassifiedLine) <- gsub("domainColumn", domainColumn, colnames(notClassifiedLine))

    abbreviations <- rbind(notClassifiedLine, abbreviations[,colnames(notClassifiedLine)])

  }



  predictionsMMFiltered <- predictionsMM %>% dplyr::filter(probability1 > probabilityThreshold)

  tumorConfusionMatrix <- caret::confusionMatrix(factor(predictionsMMFiltered$predict,
                                                 levels = base::unique(c(predictionsMM$originalCall, predictionsMM$predict))),
                                          factor(predictionsMMFiltered$originalCall, levels = base::unique(c(predictionsMM$originalCall,
                                                                                                       predictionsMM$predict))),
                                          dnn = c("Prediction", "Reference"))
  predictionFrequencies <- tumorConfusionMatrix$table %>% as.data.frame()
  predictionFrequencies$Prediction <- as.character(predictionFrequencies$Prediction)
  predictionFrequencies$Reference <- as.character(predictionFrequencies$Reference)
  missingTumors <- predictionFrequencies %>% dplyr::filter(Reference == Prediction,
                                                    Freq == 0)

  predictionFrequencies %<>% dplyr::filter(Freq != 0)

  for (j in seq(1:nrow(missingTumors))) {


    if (missingTumors$Reference[j] %notin% base::unique(predictionsMMFiltered$originalCall)) {
      missingTumors <- rbind(missingTumors, data.frame(Prediction = "Not classified",
                                                       Reference = missingTumors$Reference[j],
                                                       Freq = (nrow(predictionsMM[predictionsMM$originalCall == missingTumors$Reference[j], ]) -
                                                                 nrow(predictionsMMFiltered[predictionsMMFiltered$originalCall == missingTumors$Reference[j], ])))
      )

    }
  }

  for (j in seq(1:nrow(missingTumors))) {
    missingTumors$Prediction[j] <- abbreviations[abbreviations[,chosenClassColumn] == missingTumors$Prediction[j], "abbreviation", drop = T] %>% base::unique()
    missingTumors$Reference[j] <- abbreviations[abbreviations[,chosenClassColumn] == missingTumors$Reference[j], "abbreviation", drop = T] %>% base::unique()
  }




  difPredictions <- base::unique(predictionFrequencies$Reference) %>% as.character()


  linkClassAndDomain <- metaDataRef[ , c( chosenClassColumn, domainColumn)] %>% base::unique()

  for (i in seq(1:length(difPredictions))) {
    total <- predictionFrequencies %>% dplyr::filter(Reference == difPredictions[i])
    total$Domain <- linkClassAndDomain[linkClassAndDomain[,chosenClassColumn] == total$Reference[1], domainColumn]
    totalNum <- sum(total$Freq)
    metaTotal <-predictionsMM %>%
      dplyr::filter(originalCall == difPredictions[i]) %>%
      nrow()
    notClassified <- metaTotal - totalNum


    for (j in seq(1:nrow(total))) {
      total$Prediction[j] <- abbreviations[abbreviations[,chosenClassColumn] == total$Prediction[j], "abbreviation", drop = T] %>% unique()
      total$Reference[j] <- abbreviations[abbreviations[,chosenClassColumn] == total$Reference[j], "abbreviation", drop = T] %>% unique()
    }

    newLine <- data.frame(Prediction = "Not classified",
                          Reference = as.character(total$Reference[1]),
                          Freq = notClassified,
                          Domain = unique(total$Domain)
    )

    total <- rbind(total, newLine)

    if (i == 1) {
      confusionPlotDF <- total
    } else {
      confusionPlotDF <- rbind(confusionPlotDF, total)
    }
  }

  linkClassAndDomain$abbreviation <- NA

  for (i in seq(1:nrow(linkClassAndDomain))) {
    if (linkClassAndDomain[i, chosenClassColumn] %in% abbreviations[,chosenClassColumn]) {
      linkClassAndDomain[i, "abbreviation"] <-
        abbreviations[abbreviations[,chosenClassColumn] == linkClassAndDomain[i,chosenClassColumn], "abbreviation"] %>%
        base::unique()
    }
  }


  confusionPlotDF$Prediction <- factor(confusionPlotDF$Prediction, levels = unique(abbreviations$abbreviation[abbreviations$abbreviation %in%
                                                                                           unique(c(confusionPlotDF$Reference, confusionPlotDF$Prediction,
                                                                                                    missingTumors$Reference))]))

  missingNotClassified <- levels(confusionPlotDF$Prediction)[levels(confusionPlotDF$Prediction) %notin% confusionPlotDF[confusionPlotDF$Prediction == "Not classified","Reference"] ]
  missingNotClassified <- missingNotClassified[missingNotClassified != "Not classified"]

  for (i in seq(1:base::length(missingNotClassified))) {
    missingNotClassifiedDF <- data.frame(Prediction = "Not classified",
                                         Reference = missingNotClassified[i],
                                         Freq = 0,
                                         Domain = linkClassAndDomain[linkClassAndDomain[,"abbreviation"] == missingNotClassified[i],domainColumn]    )
    if (i == 1)  {
      missingClassifiedDF <- missingNotClassifiedDF
    } else {
      missingClassifiedDF <- rbind(missingClassifiedDF,
                                   missingNotClassifiedDF)
    }
  }

  if (subtype == T) {
    confusionPlotDF$Reference <- factor(confusionPlotDF$Reference,
                                 levels = levels(confusionPlotDF$Prediction))
  } else {
    confusionPlotDF$Reference <- factor(confusionPlotDF$Reference,
                                 levels = levels(confusionPlotDF$Prediction)[-1])
  }

  missingTumors$Domain <- NA
  for (i in seq(1:nrow(missingTumors))) {
    missingTumors$Domain[i] = linkClassAndDomain[linkClassAndDomain$abbreviation == missingTumors$Reference[i], domainColumn]
  }
  missingTumors$Prediction <- factor(missingTumors$Prediction, levels = levels(confusionPlotDF$Prediction))
  missingTumors$Reference <- factor(missingTumors$Reference, levels = levels(confusionPlotDF$Reference))

  confusionPlotDF <- rbind(confusionPlotDF, missingTumors)


    nonAvailableTiles <- getNonAvailableTiles(predictionsMM = predictionsMM,
                                              abbreviations = abbreviations,
                                              classColumn = chosenClassColumn,
                                              probabilityThreshold = probabilityThreshold)

  confusionPlotInfo <- list(confusionPlotDF = confusionPlotDF,
                            nonAvailableTiles = nonAvailableTiles,
                            abbreviations = abbreviations,
                            classColumn = classColumn,
                            higherClassColumn = higherClassColumn,
                            subtype = subtype
                            )


  return(confusionPlotInfo)
}
