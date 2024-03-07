plotTop3Accuracy <- function(meanAndSDPlotTrain,
                             meanAndSDPlotTest,
                             metaDataRef,
                             subtype
                             ) {

  meanAndSDPlotTrain$trainOrTest <- "Reference cohort"
  meanAndSDPlotTest$trainOrTest <- "Test set"

  top23 <- rbind(meanAndSDPlotTrain, meanAndSDPlotTest)
  top23$trainOrTest <- factor(top23$trainOrTest, levels = c("Reference cohort", "Test set"))

  top23Longer <- top23 %>% pivot_longer(cols = c("meanFractionCorrect",
                                                 "meanFractionCorrect2",
                                                 "meanFractionCorrect3"
  ),
  names_to = "topN",
  values_to = "fractionsCorrect"
  )
  top23Longer$fractionCorrectPercent <- paste0(round(top23Longer$fractionsCorrect * 100, 1), "%")

  top23ForSD <- top23Longer %>% pivot_longer(cols = c("sdFractionCorrect", "sdFractionCorrect2", "sdFractionCorrect3"),
                                             names_to = "whichFraction",
                                             values_to = "standardDeviation")

  #top23ForSD
  #top23ForSD$topN <- gsub("mean", "", top23ForSD$topN)
  top23ForSD$whichFraction <- gsub("sd", "mean", top23ForSD$whichFraction)
  top23ForSD %<>% filter(topN == whichFraction)

  top3Plot <- ggplot(top23Longer,
         aes(x = nCases,
             y = fractionsCorrect,
             fill = topN,
             label = fractionCorrectPercent,
             group = topN)) +
    theme_classic() +
    geom_errorbar(data = top23ForSD,
      aes(ymin= fractionsCorrect - standardDeviation, ymax= fractionsCorrect + standardDeviation), width=.2,
                  position=position_dodge(.9),
                  col = "#7F7384"

    ) +
    geom_col(position = "dodge",
             col = "black") +
    facet_grid(trainOrTest~nCases, space = "free", scales = "free" ) +
    theme(axis.text.x = element_blank(),
          legend.title = element_blank(),
          plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm")
          ) +
    geom_text(aes(y = fractionsCorrect - 0.05),
              position = position_dodge(0.9),
              size = 1.8,
              color = "white") +

    labs(y = "Fraction correct",
         x = "Number of patients per tumor type (n)"
         ) +
    scale_fill_discrete(name = "topN", labels = c("Top 1 classification",
                                                  "Top 2 classifications",
                                                  "Top 3 classifications")) +

    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

  if (subtype == T) {
    top3Plot <- top3Plot + labs(x = "Number of patients per tumor subtype (n)")

  }
  return(top3Plot)
}
