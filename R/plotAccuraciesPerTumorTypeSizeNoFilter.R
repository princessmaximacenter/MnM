#' Plot the accuracies of the different tumor type predictions per size block
#'
#' @param fractionsCorrect Dataframe containing the results per tumor type size block.
#'
#' @return ggplot object that contains the barplot with the accuracies for the
#' tumor types of different sample bins (e.g. only 3 available samples per tumor type is a bin).
#' No filtering was applied.
#' @export
#' @import ggplot2
#'
plotAccuraciesPerTumorTypeSizeNoFilter <- function(fractionsCorrect) {


  fractionsCorrect$fractionIncorrectUnfiltered <- 1 - fractionsCorrect$fractionCorrect

  figureDF <- pivot_longer(fractionsCorrect, cols = c(fractionCorrect, fractionIncorrectUnfiltered),
                            names_to = "fraction_type", values_to = "all_fractions")


  fractionsCorrect$fractionCorrectPercentUnfiltered <-  paste0("(",fractionsCorrect$fractionCorrect * 100, "%)")

  figureDF$fraction_type <- factor(figureDF$fraction_type,
                                    levels = c("fractionIncorrectUnfiltered", "fractionCorrect"))
  # Create plot
  myPlot <- ggplot(figureDF,
         aes(x = nCases,
             y = all_fractions
         )) +
    theme_classic() +
    ylab("Fraction correct and incorrect samples") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.05, hjust=1)) +
    geom_bar(aes(fill = fraction_type), stat = "identity",
             col = "grey",
             position = "stack") +

    geom_text(data = fractionsCorrect, aes(label = fractionCorrectPercentUnfiltered,
                                           y =  fractionCorrect - 0.05), size = 3) +
    geom_text(data = fractionsCorrect, aes(label = numbersClassifiedTotal, y = 0.05), size = 3.7) +
    #  geom_text(data = fractionClassified, aes(label = notClassifiedName, y = 0.98), size = 3.7) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("fractionCorrect" = "#94D6B4",
                                 "fractionIncorrectUnfiltered" =  "#F7C19C"),
                      labels = c('Incorrectly Classified', 'Correctly Classified')) +

    labs(x = "Number of patients per tumor type (n)") +
    theme(legend.title = element_blank()) +
    theme(#axis.text.x = element_text(hjust = -0.2),
      axis.title.x = element_text(vjust = -1.8),
      axis.title.y = element_text(vjust = 2)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0)))

  return(myPlot)
}
