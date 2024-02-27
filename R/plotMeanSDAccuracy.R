#' Plot the results per frequency block
#'
#' Note: This function has been replaced with the more
#'
#' @param meanNumbers
#'
#' @return
#' @export
#'
#' @examples
plotMeanSDAccuracy <- function(meanNumbers) {
meanNumbers$fractionCorrectPercent <- paste0("(",round(meanNumbers$meanFractionCorrect,2) * 100, "%)")

figureDF2 <- pivot_longer(meanNumbers, cols = c(meanFractionCorrect, meanFractionIncorrect),
                          names_to = "fraction_type", values_to = "all_fractions")

figureDF2$fraction_type <- factor(figureDF2$fraction_type,
                                  levels = c("meanFractionIncorrect", "meanFractionCorrect"))

#theme_set(theme_classic())
figureDF3 <- figureDF2 %>% filter(fraction_type == "meanFractionCorrect")
# Create plot
ggplot(figureDF2,
       aes(x = nCases,
           y = all_fractions
       )) +
  theme_classic() +
  ylab("Fraction correct and incorrect samples") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.05, hjust=1)) +
  geom_bar(aes(fill = fraction_type), stat = "identity",
           col = "darkgrey",
           position = "stack") +
  #geom_linerange(aes(ymin= minFractionCorrectMajority, ymax= maxFractionCorrectMajority), width=.2,
  #               position=position_dodge(.9),
  #              col = "darkgrey") +
  geom_errorbar(data = figureDF3, aes(ymin= all_fractions - sdFractionCorrect, ymax= all_fractions + sdFractionCorrect), width=.2,
                position=position_dodge(.9),
                col = "#7F7384") +
  geom_text(data = meanNumbers, aes(label = fractionCorrectPercent,
                                    y =  meanFractionCorrect - 0.05), size = 3) +
  geom_text(data = meanNumbers, aes(label = paste0("N = ", nSamples), y = 0.05), size = 3.7) +
  #  geom_text(data = fractionClassified, aes(label = notClassifiedName, y = 0.98), size = 3.7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("meanFractionCorrect" = "#94D6B4",
                               "meanFractionIncorrect" =  "#F7C19C"),
                    labels = c('Incorrectly Classified', 'Correctly Classified')) +

  labs(x = "Number of patients per tumor type (n)") +
  theme(legend.title = element_blank()) +
  theme(#axis.text.x = element_text(hjust = -0.2),
    axis.title.x = element_text(vjust = -1.8),
    axis.title.y = element_text(vjust = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0))) +
  geom_bar(data = figureDF3,
           aes(y = 1),
           fill = NA,
           color = "black",
           linewidth = 0.7,
           stat = "identity"
  )
}
