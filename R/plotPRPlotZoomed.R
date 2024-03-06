#' Plot precision-recall curves of classifier comparison - zoom
#'
#' @param dataPR Dataframe containing the precision and recall values at specific cutoffs,
#' for both the train and test set. Also, a minimum and maximum precision are supplied
#' to generate a shadow of the fluctuations.
#' @param otherClassifierName What is the name of the other classifier within dataPR?
#' @param removeLegend Do you want to remove the legend?
#'
#' @return Precision-recall curve for train and test dataset of M&M and other classifier, with a shadow for cross-validation.
#' This function shows a close-up of the high precision range.
plotPRPlotZoomed <- function(dataPR,
                             otherClassifierName,
                             removeLegend = T) {
dataPR$type <- factor(dataPR$type, levels = c("M&M", otherClassifierName))
dataPRFiltered <- dataPR %>% filter(trainOrTest == "Train",
                                    type == "M&M")
  PRPlotZoomed <- ggplot(dataPR,
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

    scale_color_manual(values = c("#606ca5","#a4b3f0")) +
    theme_classic() +
    theme(
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),
          axis.title = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),
                             "cm")) +

    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                       limits = c(0.85,1)
    ) +  geom_ribbon(data = dataPRFiltered,

                     aes(ymin=minPrecision, ymax=maxPrecision), alpha=0.2, colour = NA)

if (removeLegend == T) {
  PRPlotZoomed <- PRPlotZoomed + theme(legend.position = "none")
}
return(PRPlotZoomed)
}
