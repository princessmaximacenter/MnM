#' Plot UMAP tumor types
#'
#' Function to plot the UMAP coordinates for the different tumor type samples
#' belonging to one domain. All tumor types are colored differently to
#' distinguish them from one another.
#'
#' @param dataUMAP Dataframe containing the UMAP coordinates, the tumor type ($subclass),
#' the domain ($Domain) and the desired abbreviation for the tumor type ($abbreviation)
#' @param domain Which domain do you want to plot the tumor types for?
#' @param domainCol Which colors do you want to use for the tumor types?
#' @param abbreviationTumorType Abbreviations to be used in the plot for the tumor types.
#' @param useLabels Do you want to supply labels within the plot?
#' @param useManualColors Do you want to supply the colors to be used within the UMAP for the labels?
#' @return ggplot with only the datapoints of the selected domain,
#' color coded by the tumor type.
#' @import ggrepel
#' @export
#'
plotCohortTumorType <- function(dataUMAP,
                                domain,
                                domainCol,
                                classColumn,
                                abbreviationTumorType,
                                useLabels = T,
                                useManualColors = F
                                ) {



  umapDomain <- dataUMAP %>% filter(Domain == domain)
  dataLogUMAPlabels <- umapDomain %>% filter(!(duplicated(abbreviation)))

  dataLogUMAPlabels %<>% arrange(subclass)
  dataLogUMAPlabels$abbreviation <- factor(dataLogUMAPlabels$abbreviation, levels = unique(dataLogUMAPlabels$abbreviation))


  abbreviationTumorType %<>% filter(!!sym(classColumn) %in% umapDomain$subclass)
  domainTumorTypes <- abbreviationTumorType %>% filter(Domain == domain) %>% select(abbreviation) %>% deframe()


  dataLogUMAPlabels <- umapDomain %>% filter(!(duplicated(abbreviation)))

  #dataLogUMAPlabels %<>% arrange(subclass)
  #dataLogUMAPlabels$abbreviation <- factor(dataLogUMAPlabels$abbreviation, levels = unique(dataLogUMAPlabels$abbreviation))


  umapCohortTumorType <- umapDomain %>%
    ggplot(aes(x = UMAP1,
               y = UMAP2
    )) +
    theme_classic() +
    labs(color = "Tumor Type") +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =0.7),
          legend.key = element_blank()) +
    geom_point(aes(color = abbreviation),
               shape = 19,
               key_glyph = draw_key_point) +
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          axis.title = element_text(size = 14),
          axis.title.x = element_text(vjust = -1.8),
          axis.title.y = element_text(vjust = 2),
          legend.position = "none",
          plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"))

  if (useLabels == T) {
    umapCohortTumorType <- umapCohortTumorType +
      geom_label_repel(data = dataLogUMAPlabels,
                       aes(color = abbreviation,
                           label = abbreviation),
                       max.overlaps = 40,
                       size=4,
                       seed = 1,
                       label.size = 1,
                       show.legend = F,
                       fill = NA,
                       #       nudge_x = c(0.9, 3, -0.5, -3, 0.5, rep(0, times = length(unique(bloodTsne_df2$subclass))-5)),
                       #     nudge_y = c(0, -1, -2.5, rep(0, times = length(unique(bloodTsne_df2$subclass)) - 3)),
                       segment.alpha = 0
      )
  }

  if (useManualColors == T) {
    umapCohortTumorType <- umapCohortTumorType +
      scale_color_manual(values = domainCol,
                         breaks = domainTumorTypes)


  }

  return(umapCohortTumorType)


}
