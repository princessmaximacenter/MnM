#' Function to obtain (plot for) perfomance within the different domains
#'
#' @param precisionDomainTrain Dataframe containing results for reference cohort data obtained from the function 'precisionsForDomains'
#' @param precisionDomainTest Dataframe containing results for test set data obtained from the function 'precisionsForDomains'
#' @param levelsDomain What are the names you have within the domains?
#' @param returnPlot Do you want to obtain the data or (FALSE) or get the resulting plot (TRUE)?
#'
#' @return  If returnPlot = T: Plot with the performance within the different domains.
#' If returnPlot = F: Dataframe containing the mean accuracy, precision and recall within the different domains.
#' @export
#'
plotPerformanceDomain <- function(precisionDomainTrain,
         precisionDomainTest,
         levelsDomain,
         returnPlot
) {

  precisionDomainTrain$trainOrTest <- "Train"
  precisionDomainTest$trainOrTest <- "Test"

  precisionDomainsTotal <- rbind(precisionDomainTrain,
                                        precisionDomainTest)
  precisionDomainsTotal$Domain <- factor(precisionDomainsTotal$Domain, levels = levelsDomain)

  precisionDomainsTotal$trainOrTest <- factor(precisionDomainsTotal$trainOrTest, levels = c("Train", "Test"))

  precisionDomainsTotalLonger <- precisionDomainsTotal %>% tidyr::pivot_longer(cols = c(meanPrecision, meanAccuracy, meanRecall),
                                                                                      names_to = "measurementType")


  precisionDomainsTotalLonger$valuePercent <- paste0(round(precisionDomainsTotalLonger$value *100, digits = 1), "%")


  precisionDomainsTotalLonger <- precisionDomainsTotalLonger %>% tidyr::pivot_longer(cols = c("sdAccuracy", "sdPrecision", "sdRecall"),
                                                                                            names_to = "whichSD",
                                                                                            values_to = "standardDeviation")

  precisionDomainsTotalLonger$measurementType <- gsub("mean", "", precisionDomainsTotalLonger$measurementType)

  precisionDomainsTotalLonger$whichSD <- gsub("sd", "", precisionDomainsTotalLonger$whichSD)


  precisionDomainsTotalLongerFiltered <- precisionDomainsTotalLonger %>% dplyr::filter(whichSD == measurementType)

  precisionPlot <- precisionDomainsTotalLongerFiltered %>%
    ggplot2::ggplot(
      aes(x = trainOrTest,
          y= value,
          label = valuePercent,
          fill = Domain,
          alpha = trainOrTest)) +
    ggplot2::geom_col(
      color = "black") +
    theme_classic() +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05))
                       #  limits = c(0,1)
    ) +
    ggplot2::geom_errorbar(
      #data = precisionDomainsTotalLongerFiltered,
      aes(ymin = value - standardDeviation,
          ymax = value + standardDeviation),
      width=.2,
      position=position_dodge(.9),
      col = "black") +
    ggplot2::facet_grid(measurementType~Domain) +
    ggplot2::geom_text(position = position_dodge(0.9),
              aes(y = 0.5),
              size = 4,
              color = "white",
              alpha= 1
    ) +
    ggplot2::scale_alpha_manual(values = c("Train" = 1,
                                  "Test" = 0.7))+
    ggplot2::theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none",
          strip.text = element_text(
            size = 15),
          plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm")

    ) +
    ggplot2::geom_text(data = (precisionDomainsTotalLongerFiltered %>% dplyr::filter(measurementType == "Accuracy",
                                                           whichSD == "Accuracy") ),
              aes(label = paste0("N = ", nSamples), y = 1),
              size = 3.5,
              color = "black",
              alpha = 1
    )
  if (returnPlot == F) {
    return(precisionDomainsTotalLongerFiltered)
  } else {


    return(precisionPlot)

  }

}
