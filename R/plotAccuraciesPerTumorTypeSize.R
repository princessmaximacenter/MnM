#' Plot the accuracies of the filtered different tumor type predictions per size block
#'
#' @param fractionsCorrect Dataframe containing the results per tumor type size block.
#'
#' @return ggplot object that contains the barplot with the accuracies for the
#' tumor types of different sample bins (e.g. only 3 available samples per tumor type is a bin).
#'
#' @export
#' @import tidyr ggplot2
#'
plotAccuraciesPerTumorTypeSize <- function(fractionsCorrect) {

  figureDF <- tidyr::pivot_longer(fractionsCorrect, cols = c(fractionClassifiedCorrect, fractionClassifiedIncorrect, notClassified),
                           names_to = "fraction_type", values_to = "all_fractions")

  figureDF$fraction_type <- factor(figureDF$fraction_type,
                                   levels = c("notClassified","fractionClassifiedIncorrect", "fractionClassifiedCorrect"))

  fractionsCorrect$notClassifiedName <- paste0("N = ",
                                               fractionsCorrect$nSamples - fractionsCorrect$nSamplesFiltered)

  fractionsCorrect$notClassifiedNumbers <- fractionsCorrect$nSamples - fractionsCorrect$nSamplesFiltered

  fractionsCorrect2 <- tidyr::pivot_longer(fractionsCorrect, cols = c(fraction, notClassified),
                                    names_to = "fraction_type", values_to = "fraction")


  figureDF$classified <- "Yes"
  figureDF[figureDF$fraction_type == "notClassified", "classified"] <- "No"

  fractionsCorrect$fractionCorrectPercent <- paste0("(",fractionsCorrect$fractionCorrectFiltered * 100, "%)")
  fractionsCorrect$fractionCorrectPercentUnfiltered <- paste0("(", fractionsCorrect$fractionCorrect * 100, "%)")
  # Create grob
  #label <- "N = Number of\n samples per bin"
  #textgrob <- textGrob(label, gp = gpar(cex = .75), )
  #width <- unit(1, "grobwidth",textgrob) + unit(10, "points")
  #height <- unit(1, "grobheight", textgrob)+ unit(10, "points")
  #rectgrob <- rectGrob(gp=gpar(colour = "black", fill = NA), height = height, width = width)
  #labelGrob <- gTree("labelGrob", children = gList(rectgrob, textgrob))

  # Create plot
  myPlot <-  ggplot2::ggplot(figureDF,
                    aes(x = nCases,
                        y = all_fractions
                    )) +
    ggplot2::theme_classic() +
    ylab("Fraction classified samples") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.05, hjust=1)) +
    ggpattern::geom_bar_pattern(aes(fill = fraction_type, pattern = classified), stat = "identity",
                     col = "grey", position = "stack",
                     pattern_fill = "lightgrey",
                     pattern_color = "lightgrey",
                     pattern_angle = 45,
                     pattern_density = 0.025,
                     pattern_spacing = 0.03,
                     pattern_key_scale_factor = 0.6
    ) +
    ggplot2::geom_bar(data = fractionsCorrect2, aes( y = fraction),
             fill = NA,
             color = "black",
             linewidth = 0.7,
             stat = "identity"
    ) +
    ggplot2::geom_text(data = fractionsCorrect, aes(label = fractionCorrectPercent,
                                           y =  fractionClassifiedCorrect - 0.05), size = 3) +
    #geom_text(data = fractionsCorrect, aes(label = numbersClassifiedFiltered, y = 0.05), size = 3.7) +
    ggplot2::geom_text(data = fractionsCorrect, aes(label = numbersClassifiedTotal, y = 0.05), size = 3.7) +
    #geom_text(data = fractionsCorrect, aes(label = notClassifiedName, y = 0.95), size = 3.7, color = "white") +
    xlab("") + #ggtitle("Detailed diagnosis") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5),
                   legend.title = element_blank(),
                   axis.title.x = element_text(vjust = -1.8),
                   axis.title.y = element_text(vjust = 2),
                   plot.margin = unit(c(0.2, 0.6, 0.5, 0.5), "lines"),
                   legend.background = element_rect(colour = "black")
                   ) +
    ggplot2::scale_fill_manual(values = c("fractionClassifiedCorrect" = "#3ca94d",
                                 "fractionClassifiedIncorrect" =  "#ef7c6c",
                                 "notClassified" = "#4296bc"),
                      labels = c('Not Classified', 'Incorrectly Classified', 'Correctly Classified')) +
    ggpattern::scale_pattern_manual(values = c(No = "stripe", Yes = "none")) +
    ggplot2::guides(pattern = "none",
           fill = guide_legend(override.aes = list(pattern = c("stripe","none","none")))) +


    # B3C5E2
  #F7C19C
    #. 94D6B4
    # guides(color = guide_legend(override.aes = list(linetype = c(1, 0, 0) )
    # labs(x= "Tumor type sample size (n)") +
    ggplot2::labs(x = "Number of patients per tumor type (n)") +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    ggplot2::scale_x_discrete(expand = expansion(mult = c(0.1, 0)))

  return(myPlot)

}
