#' Plot F1-score and recall for separate tumors
#'
#' @param separateMeanF1 Dataframe obtained from function 'calculateSeparateF1' .
#'
#' @return Plot showing the F1 score and recall per tumor (sub)type together,
#' stratified based on the tumor's population frequency.
#' @export
#'
plotF1Recall <- function(separateMeanF1) {


  meanNumbersTrain2 <- separateMeanF1 %>% arrange(nCases, meanF1, meanRecall)
  meanNumbersTrain2$tumorType <- factor(meanNumbersTrain2$tumorType,
                                        levels = unique(meanNumbersTrain2$tumorType))


  meanNumbersTrain2 %<>% pivot_longer(cols = c("meanF1", "meanRecall"), values_to = "meanPerformance", names_to = "Splitting")
  #newDFAccuracyAndThresholdPercentage$value[is.na(newDFAccuracyAndThresholdPercentage$value)]  <- 0
  meanNumbersTrain2$Splitting <- gsub("mean","", meanNumbersTrain2$Splitting)


  meanNumbersTrain2 %>%
    ggplot(aes(x = tumorType,  y = meanPerformance, label = tumorType))  +
    theme_classic() +
    geom_col(position = "dodge", fill = "grey", #fill = "#B3C5E2",
             color = "black") +
    #geom_text(aes(y = fractionClassified + 0.005),  position = "dodge", angle = 90,
    #           hjust = 0.01) +
    facet_grid(Splitting~nCases, scales = "free_x", space = "free_x"
    ) +
    # theme(legend.position = "none",
    #       axis.text.x = element_text(angle = 90,
    #                                  hjust = 1.0),
    #       axis.ticks.x = element_blank(),
    #       axis.title = element_text(size = 14),
    #       plot.margin =  unit(c(0.8,0.8,0.8,0.8), "cm"),

    theme(axis.title.x = element_text(vjust = -1, size = 20),
          axis.title.y = element_text(vjust = 2, size = 20),

          axis.text.x = element_text(size = 15,  hjust=1, vjust = 0.05, angle = 90),
          axis.text.y = element_text(size = 15),
          strip.text = element_text(size = 11),
          legend.position = "none",
          panel.spacing = unit(1.3,"lines")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Number of patients per tumor type (n)",
         y = "Performance Unfiltered Data")



}
