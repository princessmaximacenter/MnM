plotSampleHeterogeneity <- function(tumorHeterogeneityTrain,
                                    tumorHeterogeneityTest,
                                    FFPEResults,
                                    FFPEResultsSubtype,
                                    subtype

                                    ) {

  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(10, "Spectral"))


  FFPEResults$type <- "TumorType"
  FFPEResults$trainOrTest <- "Test"

  FFPEResults <- FFPEResults %>% group_by(labelType,
                                                  trainOrTest,
                                                  type,
                                                  numberSamples
  ) %>%
    summarise(
      meanAccuracy = mean(accuracy),
      meanPrecision = mean(precision),
      meanRecall = mean(recall),
      #meanIncorrect = mean(1 - fractionCorrect),
      #meanFractionIncorrectFiltered = mean(1 - fractionCorrectFiltered),
      sdAccuracy = sd(accuracy),
      sdPrecision = sd(precision),
      sdRecall = sd(recall)
    )

  # FFPEResults <- FFPEResults[c("labelType", "trainOrTest",
  #                              "type", "numberSamples",
  #                              "meanAccuracy", "meanPrecision",
  #                              "meanRecall", "sdAccuracy",
  #                              "sdPrecision")]

  FFPEResultsSubtype$type <- "TumorSubtype"
  FFPEResultsSubtype$trainOrTest <- "Test"

  FFPEResultsSubtype  <- FFPEResultsSubtype %>% group_by(labelType,
                                                  trainOrTest,
                                                  type,
                                                  numberSamples
  ) %>%
    summarise(
      meanAccuracy = mean(accuracy),
      meanPrecision = mean(precision),
      meanRecall = mean(recall),
      #meanIncorrect = mean(1 - fractionCorrect),
      #meanFractionIncorrectFiltered = mean(1 - fractionCorrectFiltered),
      sdAccuracy = sd(accuracy),
      sdPrecision = sd(precision),
      sdRecall = sd(recall)
    )
  totalHeterogeneityDF <- rbind(tumorHeterogeneityTrain,
                                tumorHeterogeneityTest,
                                FFPEResults,
                                FFPEResultsSubtype
                                )

  totalHeterogeneityDF$labelType <- as.character(totalHeterogeneityDF$labelType)

  totalHeterogeneityDF$labelType <- factor(totalHeterogeneityDF$labelType,
                                            levels = c("All", "Primary",
                                                       "Recurrence",
                                                       "Metastasis",
                                                       "No treatment",
                                                       "Systemic treatment",
                                                       "FFPE"))

  colorBrewerColors <-  getPalette(length(unique(totalHeterogeneityDF$labelType)))

  tumorHeterogeneityLonger <- totalHeterogeneityDF %>% pivot_longer(cols = c("meanAccuracy", "meanPrecision", "meanRecall"),
                                                                    names_to = "measurementType")
  tumorHeterogeneityLonger$valuePercent <- paste0(round(tumorHeterogeneityLonger$value * 100,1), "%")

  tumorHeterogeneityLonger2 <- tumorHeterogeneityLonger %>% pivot_longer(cols = c("sdAccuracy", "sdPrecision", "sdRecall"),
                                                                         names_to = "whichSD",
                                                                         values_to = "standardDeviation")

  tumorHeterogeneityLonger2$measurementType <- gsub("mean", "", tumorHeterogeneityLonger2$measurementType)
 # tumorHeterogeneityLonger2$measurementType <- gsub("meanPrecision", "precision", tumorHeterogeneityLonger2$measurementType)
 # tumorHeterogeneityLonger2$measurementType <- gsub("meanRecall", "recall", tumorHeterogeneityLonger2$measurementType)

 # tumorHeterogeneityLonger2$whichSD <- gsub("sdAccuracy", "Accuracy", tumorHeterogeneityLonger2$whichSD)
 # tumorHeterogeneityLonger2$whichSD <- gsub("sdPrecision", "recision", tumorHeterogeneityLonger2$whichSD)
  tumorHeterogeneityLonger2$whichSD <- gsub("sd", "", tumorHeterogeneityLonger2$whichSD)

  tumorHeterogeneityLonger3 <- tumorHeterogeneityLonger2 %>% filter(whichSD == measurementType)

  tumorHeterogeneityLonger3$trainOrTest <-factor(tumorHeterogeneityLonger3$trainOrTest,
                                                 levels = c("Train", "Test"))

  if (subtype == F) {
  tumorHeterogeneityLonger3  %<>% filter(type == "TumorType")
    } else {
      tumorHeterogeneityLonger3  %<>% filter(type == "TumorSubtype")
    }
      tumorHeterogeneityLonger3    %>%
    ggplot(
      aes(x = trainOrTest,
          y = value,
          fill = labelType,
          #color = TrainOrTest,

          # group = type

      )) +
    theme_classic() +
    geom_errorbar(aes(ymin = value - standardDeviation,
                      ymax = value + standardDeviation),
                  width=.7,
                  position=position_dodge(.9),
                  col = "#7F7384") +
    geom_col(
      position = position_dodge(0.9),
      #position = "dodge",
      # stat = "identity"
      color = "black"
    ) +
    # scale_fill_manual(values = c("Train" = "#606ca5",
    #                              "Test" =  "#a4b3f0")) +
    #scale_color_manual(values = c("train" = "black",
    #                               "test" = "black")) +
    facet_grid(measurementType  ~ labelType) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       limits = c(0,1.10)
    ) +
    geom_text(position = position_dodge(0.9),
              aes(y = 0.5, label = valuePercent),
              size = 2.5,
              color = "white",

    ) +
        scale_fill_manual(values=colorBrewerColors) +


    geom_text(data = (tumorHeterogeneityLonger3 %>% filter(measurementType == "Accuracy")),
              aes(label = paste0("N = ", numberSamples), y = value + 0.05),
              size = 2.5,
              color = "black") +
    theme(axis.title.y = element_blank() ,
          axis.title.x = element_blank(),
          legend.position = "none"
    )

}
