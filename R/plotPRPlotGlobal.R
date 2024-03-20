#' Plot precision-recall curves of classifier comparison - global
#'
#' @param dataPR Dataframe containing the precision and recall values at specific cutoffs,
#' for both the train and test set. Also, a minimum and maximum precision are supplied
#' to generate a shadow of the fluctuations.
#' @param otherClassifierName What is the name of the other classifier within dataPR?
#' @param removeLegend Do you want to remove the legend?
#'
#' @return Precision-recall curve for train and test dataset of M&M and other classifier, with a shadow for cross-validation.
#' This function shows the full range from 0-1 for the precision values.
#'
plotPRPlotGlobal <- function(dataPR,
                             otherClassifierName) {
  dataPR$type <- factor(dataPR$type, levels = c("M&M", otherClassifierName))

  dataPRFiltered <- dataPR %>% filter(trainOrTest == "Train",
                                      type == "M&M")
  ggplot(dataPR,
         aes(x = Recall,
             linetype = trainOrTest,
             color = type,
             y = Precision
         )) +
    scale_x_reverse(limits = c(1,0), expand = expansion(mult = c(0,0)),
                    labels = scales::number_format(accuracy = 0.01)) +
    geom_line() +
    # scale_color_manual(values = c("M&M" = "#606ca5",
    #                               "Other" =  "#a4b3f0")) +

    scale_color_manual(values = c("#606ca5","#bdc6e5")) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),
                             "cm")) +

    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                       limits = c(0.0,1)
    ) +
    geom_ribbon(data = dataPRFiltered,

                aes(ymin=minPrecision, ymax=maxPrecision), alpha=0.2, colour = NA)

}
