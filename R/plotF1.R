#' Plot F1-scores for tumors from different population frequencies
#'
#' @param dataMeanF1 Dataframe obtained from function 'calculateMeanAndSDAccuracy'.
#' @param separateMeanF1 Dataframe obtained from function 'calculateSeparateF1' .
#' @param subtype Do you want to visualize the tumor subtype (TRUE) or tumor type (FALSE)?
#' dataMeanF1 and separateMeanF1 should have been generated on the same level,
#' being either subtype or tumor type.
#'
#' @return Plot showing the average F1 values for tumors stratified by population frequency,
#' and separate F1 values for the separate tumors within the population frequency blocks.
#' @export
plotF1 <- function(dataMeanF1,
                   separateMeanF1,
                   subtype,
                   whichSeed =1
                   ) {

  set.seed(whichSeed)
  separateMeanF1 %<>% filter(!is.na(meanF1))

  dataMeanF1$meanF1Percent <- paste0(round(dataMeanF1$meanF1,3) * 100, "%")
  dataMeanF1$sdF1[is.na(dataMeanF1$sdF1)] <- 0
  ourPlot <- dataMeanF1 %>%


    ggplot(
      aes(x = type,
          y = meanF1
      )) +
    theme_classic() +
    ylab("Mean F1-score") +
    geom_rect(aes(
      xmin = 0,
      xmax = 2,
      ymin = -Inf,
      ymax = Inf,
      fill = nCases
    ),
    alpha = .2
    ) +
    geom_hpline(#aes(fill = fraction_type),
      width = 2,
      size = 1,
      stat = "identity",
      col = "darkgrey",
      #col = "black",
      # fill = "grey",
      position = "stack") +
    geom_text(#data = dataMeanF1,
              aes(label = meanF1Percent,
                                     y =  meanF1 - sdF1 - 0.04), size = 5) +

    geom_jitter(data = separateMeanF1, aes(y = meanF1),
                shape = 19, alpha = 0.7) +
    geom_errorbar(#data = figureDF3,
      aes(ymin= meanF1 - sdF1, ymax= meanF1 + sdF1), width=.2,
      position=position_dodge(.9),
      #col = "#7F7384"
      col = "black"

    ) +
    scale_fill_manual(
      values = c(
        '#7a5d7e',
        'lightgrey',
        '#7a5d7e',
        'lightgrey',
        '#7a5d7e',
        'lightgrey',
        '#7a5d7e',
        'lightgrey'
      )) +
    facet_grid(~nCases, scales = "free_x", space = "free_x") +
    labs(x = "Number of patients per tumor type (n)") +

    theme(#axis.text.x = element_text(hjust = -0.2),
      legend.title = element_blank(),
      panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
          axis.title.x = element_text(vjust = -1, size = 20),
          axis.title.y = element_text(vjust = 2, size = 20),

          axis.text.x = element_text(size = 15,  hjust=1, vjust = 0.05, angle = 90),
          axis.text.y = element_text(size = 15),
          strip.text = element_text(size = 11),
          legend.position = "none",
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +

    scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0)))


  if (subtype == T) {
    ourPlot <- ourPlot + labs(x = "Number of patients per tumor subtype (n)")
  }
  return(ourPlot)
}
