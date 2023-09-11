#' Title
#'
#' @param predictionsMM
#' @param metaDataRef
#' @param abbreviations
#' @param classColumn
#' @param defineTumorWithColor
#' @param probabilityScore
#' @import tidyverse dplyr magrittr
#' @return
#' @export
#'
#' @examples
getConfusionMatrixPlot <- function(predictionsMM,
         metaDataRef,
         abbreviations,
         classColumn,
         defineTumorWithColor = F,
         probabilityScore = 0.8) {
  predictionsMMFiltered <- predictionsMM %>% filter(probability1 > probabilityScore)
  tumorConfusionMatrix <- confusionMatrix(factor(predictionsMMFiltered$predict,
                                                 levels = unique(c(predictionsMMFiltered$originalCall, predictionsMMFiltered$predict))),
                                          factor(predictionsMMFiltered$originalCall, levels = unique(c(predictionsMMFiltered$originalCall,
                                                                                                       predictionsMMFiltered$predict))),
                                          dnn = c("Prediction", "Reference"))
  predictionFrequencies <- tumorConfusionMatrix$table %>% as.data.frame() #%>% filter(Freq != 0)
  predictionFrequencies$Prediction <- as.character(predictionFrequencies$Prediction)
  predictionFrequencies$Reference <- as.character(predictionFrequencies$Reference)
  missingTumors <- predictionFrequencies %>% filter(Reference == Prediction,
                                                    Freq == 0)

  for (j in seq(1:nrow(missingTumors))) {
    missingTumors$Prediction[j] <- abbreviations[abbreviations[,classColumn] == missingTumors$Prediction[j], "abbreviation", drop = T]
    missingTumors$Reference[j] <- abbreviations[abbreviations[,classColumn] == missingTumors$Reference[j], "abbreviation", drop = T]
  }


  predictionFrequencies %<>% filter(Freq != 0)

  difPredictions <- unique(predictionFrequencies$Reference) %>% as.character()


  linkClassAndDomain <- metaDataRef[ , c( classColumn, "Domain")] %>% unique

  for (i in seq(1:length(difPredictions))) {
    total <- predictionFrequencies %>% filter(Reference == difPredictions[i])
    total$Domain <- linkClassAndDomain[linkClassAndDomain[,classColumn] == total$Reference[1], "Domain"]
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
      ggplotDF <- total
    } else {
      ggplotDF <- rbind(ggplotDF, total)
    }
  }

  linkClassAndDomain$abbreviation <- NA

  for (i in seq(1:nrow(linkClassAndDomain))) {
    if (linkClassAndDomain[i, classColumn] %in% abbreviations[,classColumn]) {
      linkClassAndDomain[i, "abbreviation"] <- abbreviations[abbreviations[,classColumn] == linkClassAndDomain[i,classColumn], "abbreviation"]
    }
  }


  ggplotDF$Prediction <- factor(ggplotDF$Prediction, levels = abbreviations$abbreviation[abbreviations$abbreviation %in%
                                                                                           unique(c(ggplotDF$Reference, ggplotDF$Prediction))])

  missingNotClassified <- levels(ggplotDF$Prediction)[levels(ggplotDF$Prediction) %notin% ggplotDF[ggplotDF$Prediction == "Not classified","Reference"] ]
  missingNotClassified <- missingNotClassified[missingNotClassified != "Not classified"]

  for (i in seq(1:length(missingNotClassified))) {
    missingNotClassifiedDF <- data.frame(Prediction = "Not classified",
                                         Reference = missingNotClassified[i],
                                         Freq = 0,
                                         Domain = linkClassAndDomain[linkClassAndDomain[,"abbreviation"] == missingNotClassified[i],"Domain"]    )
    if (i == 1)  {
      missingClassifiedDF <- missingNotClassifiedDF
    } else {
      missingClassifiedDF <- rbind(missingClassifiedDF,
                                   missingNotClassifiedDF)
    }
  }
  ggplotDF <- rbind(ggplotDF, missingClassifiedDF)

  if (defineTumorWithColor == T) {
    ggplotDF$Reference <- factor(ggplotDF$Reference,
                                 levels = levels(ggplotDF$Prediction))
  } else {
    ggplotDF$Reference <- factor(ggplotDF$Reference,
                                 levels = levels(ggplotDF$Prediction)[-1])
  }

  missingTumors$Domain <- NA
  for (i in seq(1:nrow(missingTumors))) {
    missingTumors$Domain[i] = linkClassAndDomain[linkClassAndDomain$abbreviation == missingTumors$Prediction[i], "Domain"]
  }
  missingTumors$Prediction <- factor(missingTumors$Prediction, levels = levels(ggplotDF$Prediction))
  missingTumors$Reference <- factor(missingTumors$Reference, levels = levels(ggplotDF$Reference))

  ggplotDF <- rbind(ggplotDF, missingTumors)

  return(ggplotDF)
}
