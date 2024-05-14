#' Plot F1-scores and recall per tumor frequency
#'
#' The plotF1Recall function can be used in two ways, to generate different plot types.
#' One plot shows the tumor population group frequency averages for F1 scores, together with the individual F1-values.
#' No recall is displayed in this plot.
#'
#' The second plot shows the overall F1-scores for both confident and non-confident classifications together with the associated recall rates.
#' This type of plot can be created by only specifying 'separateMeanF1'.
#'
#' @param separateMeanF1 Dataframe obtained from function 'calculateSeparateF1'.
#' Contains all the calculated F1 scores and associated recall rates per tumor (sub)type.
#' For this setup, filterOrNot within 'calculateSeparateF1' should be put to FALSE.
#' @param dataMeanF1 Dataframe obtained from function 'calculateMeanAndSDAccuracy'.
#' Contains all the calculated F1 scores and associated recall rates stratified by population frequency.
#' @param whichSeed For reproducibility, the seed can be specified with this parameter.
#' @return If dataMeanF1 is not specified: Plot showing the F1 score and recall per tumor (sub)type together,
#' stratified based on the tumor's population frequency.
#'
#' If dataMeanF1 is specified: Plot showing the tumor population group frequency averages for F1 scores,
#' together with the individual values of the tumor (sub)types that are part of the frequency group.
#'
#' @export
#'
plotF1Recall <- function(separateMeanF1,
                         dataMeanF1 = NA,
                         whichSeed = 1) {


  if (require("ungeviz") == F) {
    remotes::install_github("fwallis/ungeviz")
  }

  set.seed(whichSeed)
  separateMeanF1 %<>% dplyr::filter(!is.na(meanF1))
  meanNumbersTrain2 <- separateMeanF1 %>% dplyr::arrange(nCases, meanF1, meanRecall)
  meanNumbersTrain2$tumorType <- factor(meanNumbersTrain2$tumorType,
                                        levels = unique(meanNumbersTrain2$tumorType))

  subtype <- separateMeanF1$subtype[1]

  meanNumbersTrain2 %<>% tidyr::pivot_longer(cols = c("meanF1", "meanRecall"), values_to = "meanPerformance", names_to = "Splitting")

    meanNumbersTrain2$Splitting <- gsub("mean","", meanNumbersTrain2$Splitting)

if (is.na(dataMeanF1)[1]) {
  recallPlot <- meanNumbersTrain2 %>%
    ggplot2::ggplot(ggplot2::aes(x = tumorType,  y = meanPerformance, label = tumorType))  +
    ggplot2::theme_classic() +
    ggplot2::geom_col(position = "dodge", fill = "grey",
             color = "black") +
    ggplot2::facet_grid(Splitting~nCases, scales = "free_x", space = "free_x") +
    ggplot2::theme(axis.title.x = ggplot2::element_text(vjust = -1, size = 20),
          axis.title.y = ggplot2::element_text(vjust = 2, size = 20),
          axis.text.x = ggplot2::element_text(size = 15,  hjust=1, vjust = 0.05, angle = 90),
          axis.text.y = ggplot2::element_text(size = 15),
          strip.text = ggplot2::element_text(size = 11),
          legend.position = "none",
          panel.spacing = ggplot2::unit(1.3,"lines")
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = "Number of patients per tumor type (n)",
         y = "Performance Unfiltered Data")
} else {

  dataMeanF1$meanF1Percent <- base::paste0(base::round(dataMeanF1$meanF1,3) * 100, "%")
  dataMeanF1$sdF1[base::is.na(dataMeanF1$sdF1)] <- 0
  recallPlot <- dataMeanF1 %>%

    ggplot2::ggplot(
      ggplot2::aes(x = type,
          y = meanF1
      )) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Mean F1-score") +
    ggplot2::geom_rect(
      ggplot2::aes(
      xmin = 0,
      xmax = 2,
      ymin = -Inf,
      ymax = Inf,
      fill = nCases
    ),
    alpha = .2
    ) +
    ungeviz::geom_hpline(#aes(fill = fraction_type),
      width = 2,
      size = 1,
      stat = "identity",
      col = "darkgrey",
      #col = "black",
      # fill = "grey",
      position = "stack") +
    ggplot2::geom_text(#data = dataMeanF1,
      ggplot2::aes(label = meanF1Percent,
          y =  meanF1 - sdF1 - 0.04), size = 5) +

    ggplot2::geom_jitter(data = separateMeanF1,
                         ggplot2::aes(y = meanF1),
                shape = 19, alpha = 0.7) +
    ggplot2::geom_errorbar(#data = figureDF3,
      ggplot2::aes(ymin= meanF1 - sdF1, ymax= meanF1 + sdF1), width=.2,
      position=ggplot2::position_dodge(.9),
      #col = "#7F7384"
      col = "black"

    ) +
    ggplot2::scale_fill_manual(
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
    ggplot2::facet_grid(~nCases, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = "Number of patients per tumor type (n)") +

    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(vjust = -1, size = 20),
      axis.title.y = ggplot2::element_text(vjust = 2, size = 20),

      axis.text.x = ggplot2::element_text(size = 15,  hjust=1, vjust = 0.05, angle = 90),
      axis.text.y = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(size = 11),
      legend.position = "none",
      plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +

    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.1, 0)))


}
  if (subtype == T) {
    recallPlot <- recallPlot + ggplot2::labs(x = "Number of patients per tumor subtype (n)")
  }


return(recallPlot)


}
