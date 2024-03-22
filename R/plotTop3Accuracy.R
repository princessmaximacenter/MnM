#' Plot accuracy within the top 3 highest scoring classification labels M&M
#'
#' @param meanAndSDPlotTrain Dataframe containing the average performance of tumor classifications within a certain frequency range (nCases) for the reference cohort.
#' Object obtained from running the function calculateMeanAndSDAccuracy
#' @param meanAndSDPlotTest  Dataframe containing the average performance of tumor classifications within a certain frequency range (nCases) for the test set.
#' Object obtained from running the function calculateMeanAndSDAccuracy
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the reference cohort.
#' @param subtype Do you want to analyse the results on the tumor subtype classification level?
#'
#' @return Plot showing the accuracy if not only the top scoring classification label, but also the second and third highest labels are taken into consideration.
#' Plot is segregated into the different class frequencies.
#' @export
#'
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

  top23Longer <- top23Longer %>% pivot_longer(cols = c("sdFractionCorrect", "sdFractionCorrect2", "sdFractionCorrect3"),
                                             names_to = "whichFraction",
                                             values_to = "standardDeviation")

  #top23ForSD
  #top23ForSD$topN <- gsub("mean", "", top23ForSD$topN)
  top23Longer$whichFraction <- gsub("sd", "mean", top23Longer$whichFraction)
  top23Longer %<>% filter(topN == whichFraction)
  top23Longer$topN[which(top23Longer$topN == "meanFractionCorrect")] <- "Top 1"
  top23Longer$topN <- gsub("meanFractionCorrect", "Top ", top23Longer$topN)


  top3Plot <- ggplot(top23Longer,
         aes(x = topN,
             y = fractionsCorrect,
             #fill = topN,
             alpha = topN,
             label = fractionCorrectPercent,
             group = topN)) +
    theme_classic() +

    geom_errorbar(
      aes(ymin= fractionsCorrect - standardDeviation, ymax= fractionsCorrect + standardDeviation), width=.2,
      alpha = 1,
                  position=position_dodge(.9),
                  col = "black"

    ) +
    geom_col(position = "dodge",
             col = "black") +
    facet_grid(trainOrTest~nCases, space = "free", scales = "free" ) +
    theme(#axis.text.x = element_text(angle = 90),
          legend.position = "none",
          #plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),

          axis.title.x = element_text(vjust = -1, size = 15),
          axis.title.y = element_text(vjust = 2, size = 15),

          axis.text.x = element_text(size = 10,  hjust=1, vjust = 0.05, angle = 90),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 11),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +
    geom_text(aes(y = fractionsCorrect + 0.05),
              position = position_dodge(0.9),
              size = 2.5,
              alpha = 1,
              color = "black") +

    labs(y = "Fraction correct classifications",
         x = "Number of patients per tumor type (n)"
         ) +
    scale_fill_discrete(name = "topN", labels = c("Top 1 classification",
                                                  "Top 2 classifications",
                                                  "Top 3 classifications")) +

    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       limits = c(0,1.05))

  if (subtype == T) {
    top3Plot <- top3Plot + labs(x = "Number of patients per tumor subtype (n)")

  }
  return(top3Plot)
}
