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

    scale_color_manual(values = c("#606ca5","#a4b3f0")) +
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
