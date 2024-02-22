plotComparisonBetweenClassifiers <- function(comparisonDF) {

  comparisonDF <- comparisonDF %>% group_by(Subset, TrainOrTest, type) %>%
    summarise(
      numberSamples = mean(nSamples),
      meanAccuracy = mean(accuracy),
      meanPrecision = mean(precision),
      meanRecall = mean(recall),
      sdAccuracy = sd(accuracy),
      sdPrecision = sd(precision),
      sdRecall = sd(recall)
    )


  comparisonDFLonger <- pivot_longer(comparisonDF, cols = c(meanAccuracy,
                                                                      meanPrecision,
                                                                      meanRecall
  ), names_to = "measurementType",
  values_to = "value")
  #%>% select(all_of(Subset, TrainOrTest, measurementType))

  comparisonDFLonger <- comparisonDFLonger %>% pivot_longer(cols = c("sdAccuracy", "sdPrecision", "sdRecall"),
                                                                   names_to = "whichSD",
                                                                   values_to = "standardDeviation")

  comparisonDFLonger$measurementType <- gsub("mean", "", comparisonDFLonger$measurementType)

  comparisonDFLonger$whichSD <- gsub("sd", "", comparisonDFLonger$whichSD)


  comparisonDFLonger$TrainOrTest <- factor(comparisonDFLonger$TrainOrTest,
                                               levels = c("Train", "Test"))
  comparisonDFLonger$valuePercent <- paste0(round(comparisonDFLonger$value *100, digits = 1), "%")

  comparisonDFLonger <- comparisonDFLonger %>% filter(whichSD == measurementType)


  ggplot(comparisonDFLonger,
         aes(x = TrainOrTest,
             y = value,
             fill = type,
             #color = TrainOrTest,
             label = valuePercent,
             group = type)) +
    theme_classic() +
    geom_errorbar(#data = tTestDFCombinedLonger3,
      aes(ymin= value - standardDeviation, ymax= value + standardDeviation), width=.2,
                  position=position_dodge(.9),
                  col = "black"

    ) +
    geom_col(
      position = position_dodge(0.9),
      #position = "dodge",
      # stat = "identity"
      color = "black"
    ) +
    scale_fill_manual(values = c("M&M" = "#606ca5",
                                 "Other" =  "#a4b3f0")) +
    #scale_color_manual(values = c("train" = "black",
    #                               "test" = "black")) +
    facet_grid(measurementType  ~ Subset) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       limits = c(0,1.05)) +
    geom_text(position = position_dodge(0.9),
              aes(y = 0.5),
              size = 2.8,
              color = "white"
    ) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 13),
          legend.position = "none",
          strip.text = element_text(
            size = 15)

    ) +
    geom_text(data = (comparisonDFLonger %>% filter(type == "M&M",
                                                        measurementType == "Accuracy")),
              aes(label = paste0("N = ", numberSamples), y = 1.04),
              size = 3.5,
              color = "black")

}
