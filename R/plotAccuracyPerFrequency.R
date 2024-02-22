plotAcccuracyPerFrequency <- function(meanAndSDPlotTrain,
                                      meanAndSDPlotTest,
                                      subtype ) {
  meanAndSDPlotTrain$type <- "Train"

  meanAndSDPlotTest$type <- "Test"
    totalMeanThing <- rbind(meanAndSDPlotTrain,
                          meanAndSDPlotTest)

    totalMeanThing$type <- factor( totalMeanThing$type, levels = c("Train", "Test"))
  totalMeanThing$fractionCorrectPercent <- paste0("(",round(totalMeanThing$meanFractionCorrectFiltered,2) * 100, "%)")

  figureDF2 <- pivot_longer(totalMeanThing, cols = c(meanFractionCorrectFiltered, meanFractionIncorrectFiltered),
                            names_to = "fraction_type", values_to = "all_fractions")

  figureDF2$fraction_type <- factor(figureDF2$fraction_type,
                                    levels = c("meanFractionIncorrectFiltered", "meanFractionCorrectFiltered"))

  #theme_set(theme_classic())
  figureDF3 <- figureDF2 %>% filter(fraction_type == "meanFractionCorrectFiltered")
  # Create plot
  ourImage <- ggplot(figureDF2,
         aes(x = type,
             y = all_fractions
         )) +
    theme_classic() +
    ylab("Fraction correct and incorrect samples") +
   # theme(axis.text.x = element_text(angle = 90, vjust = 0.05, hjust=1)) +
    geom_bar(aes(fill = fraction_type),
             width = 1,
             stat = "identity",
             #col = "darkgrey",
             col = "black",
             position = "stack") +
    #geom_linerange(aes(ymin= minFractionCorrectMajority, ymax= maxFractionCorrectMajority), width=.2,
    #               position=position_dodge(.9),
    #              col = "darkgrey") +
    geom_errorbar(data = figureDF3, aes(ymin= all_fractions - sdFractionCorrect, ymax= all_fractions + sdFractionCorrect), width=.2,
                  position=position_dodge(.9),
                  col = "#7F7384"

    ) +
    geom_text(data = totalMeanThing, aes(label = fractionCorrectPercent,
                                         y =  meanFractionCorrectFiltered - 0.05), size = 3) +
    geom_text(data = totalMeanThing, aes(label = paste0("N = ", meanCasesFiltered), y = 0.05), size = 3.5,
              angle = 90) +
    #  geom_text(data = fractionClassified, aes(label = notClassifiedName, y = 0.98), size = 3.7) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("meanFractionCorrectFiltered" = "#94D6B4",
                                 "meanFractionIncorrectFiltered" =  "#F7C19C"),
                      labels = c('Incorrectly Classified', 'Correctly Classified')) +

    labs(x = "Number of patients per tumor type (n)") +
    theme(legend.title = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_blank(),
      axis.title.x = element_text(vjust = -1.8, size = 25),
      axis.title.y = element_text(vjust = 2, size = 25),
      axis.text = element_text(size = 13, vjust = 0.05),
                               panel.border = element_blank(), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"),
                               panel.spacing = unit(1.2,"lines"),
                               plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm")
#, hjust = 1)
     # axis.text.x = element_text(angle = 90, vjust = 0.05, hjust=1)

      ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0))) +
    # geom_bar(data = figureDF3,
    #        aes(y = 1),
    #           fill = NA,
    #           color = "black",
    #           linewidth = 0.7,
    #           width = 1,
    #           stat = "identity") +
    facet_grid(~nCases, scales = "free_x", space = "free_x")
  if (subtype == T) {
    ourImage <- ourImage + labs(x = "Number of patients per tumor subtype (n)")

  }
ourImage
}
