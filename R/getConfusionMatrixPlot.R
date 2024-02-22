#' Get confusion matrix information
#'
#' This function is designed to extract which tiles of the confusion matrix
#' do contain reference-prediction information.
#'
#' @param predictionsMM Dataframe showing the top 3 predictions for the tumor (sub)type, together with their probability scores.
#' @param metaDataRef Metadata file containing the links between the patients and
#' the tumor (sub)type diagnosis within the reference cohort.
#' @param abbreviations Dataframe containing the links between the tumor (sub)type,
#' the abbreviation required in the plot, and the domain.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param domainColumn Column in the metadata file that contains the tumor domain labels.
#' @param defineTumorWithColor Do you want to have an extra column where the tumors can
#' be specified with a color?
#' @param probabilityScore Which probability score do you want to use as a cutoff
#' for being classified or non-classified?
#' @import tidyverse dplyr magrittr
#' @return Dataframe specifying how often certain reference-prediction
#' combinations are present.
#' @export
#'
getConfusionMatrixPlot <- function(predictionsMM,
         metaDataRef,
         abbreviations,
         classColumn,
         domainColumn,
         defineTumorWithColor = F,
         probabilityScore = 0.8) {
  predictionsMMFiltered <- predictionsMM %>% filter(probability1 > probabilityScore)
  #tumorConfusionMatrix <- confusionMatrix(factor(predictionsMMFiltered$predict,
   #                                              levels = unique(c(predictionsMMFiltered$originalCall, predictionsMMFiltered$predict))),
    #                                      factor(predictionsMMFiltered$originalCall, levels = unique(c(predictionsMMFiltered$originalCall,
     #                                                                                                  predictionsMMFiltered$predict))),
      #                                    dnn = c("Prediction", "Reference"))

  tumorConfusionMatrix <- confusionMatrix(factor(predictionsMMFiltered$predict,
                                                 levels = unique(c(predictionsMM$originalCall, predictionsMM$predict))),
                                          factor(predictionsMMFiltered$originalCall, levels = unique(c(predictionsMM$originalCall,
                                                                                                       predictionsMM$predict))),
                                          dnn = c("Prediction", "Reference"))
  predictionFrequencies <- tumorConfusionMatrix$table %>% as.data.frame() #%>% filter(Freq != 0)
  predictionFrequencies$Prediction <- as.character(predictionFrequencies$Prediction)
  predictionFrequencies$Reference <- as.character(predictionFrequencies$Reference)
  missingTumors <- predictionFrequencies %>% filter(Reference == Prediction,
                                                    Freq == 0)

  predictionFrequencies %<>% filter(Freq != 0)

  for (j in seq(1:nrow(missingTumors))) {


    if (missingTumors$Reference[j] %notin% unique(predictionsMMFiltered$originalCall)) {
      missingTumors <- rbind(missingTumors, data.frame(Prediction = "Not classified",
                                                       Reference = missingTumors$Reference[j],
                                                       Freq = (nrow(predictionsMM[predictionsMM$originalCall == missingTumors$Reference[j], ]) -
                                                                 nrow(predictionsMMFiltered[predictionsMMFiltered$originalCall == missingTumors$Reference[j], ])))
      )

    }
  }


  for (j in seq(1:nrow(missingTumors))) {
    missingTumors$Prediction[j] <- abbreviations[abbreviations[,classColumn] == missingTumors$Prediction[j], "abbreviation", drop = T]
    missingTumors$Reference[j] <- abbreviations[abbreviations[,classColumn] == missingTumors$Reference[j], "abbreviation", drop = T]

  }




  difPredictions <- unique(predictionFrequencies$Reference) %>% as.character()


  linkClassAndDomain <- metaDataRef[ , c( classColumn, domainColumn)] %>% unique

  for (i in seq(1:length(difPredictions))) {
    total <- predictionFrequencies %>% filter(Reference == difPredictions[i])
    total$Domain <- linkClassAndDomain[linkClassAndDomain[,classColumn] == total$Reference[1], domainColumn]
    totalNum <- sum(total$Freq)
    metaTotal <-predictionsMM %>%
      filter(originalCall == difPredictions[i]) %>%
      nrow()
    notClassified <- metaTotal - totalNum


    for (j in seq(1:nrow(total))) {
      total$Prediction[j] <- abbreviations[abbreviations[,classColumn] == total$Prediction[j], "abbreviation", drop = T]
      total$Reference[j] <- abbreviations[abbreviations[,classColumn] == total$Reference[j], "abbreviation", drop = T]
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
    if (linkClassAndDomain[i, classColumn] %in% abbreviations[,classColumn]) {
      linkClassAndDomain[i, "abbreviation"] <- abbreviations[abbreviations[,classColumn] == linkClassAndDomain[i,classColumn], "abbreviation"]
    }
  }


  confusionPlotDF$Prediction <- factor(confusionPlotDF$Prediction, levels = abbreviations$abbreviation[abbreviations$abbreviation %in%
                                                                                           unique(c(confusionPlotDF$Reference, confusionPlotDF$Prediction,
                                                                                                    missingTumors$Reference))])

  missingNotClassified <- levels(confusionPlotDF$Prediction)[levels(confusionPlotDF$Prediction) %notin% confusionPlotDF[confusionPlotDF$Prediction == "Not classified","Reference"] ]
  missingNotClassified <- missingNotClassified[missingNotClassified != "Not classified"]

  for (i in seq(1:length(missingNotClassified))) {
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
  #confusionPlotDF <- rbind(confusionPlotDF, missingClassifiedDF)

  if (defineTumorWithColor == T) {
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

  return(confusionPlotDF)
}
