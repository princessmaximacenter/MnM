#'  Plot UMAP Domains
#'
#' Function to plot the UMAP coordinates for the different domains.
#' All domains are colored differently to
#' distinguish them from one another.
#'
#' @param dataUMAP Dataframe containing the UMAP coordinates, the tumor type ($subclass),
#' the domain ($Domain) and the desired abbreviation for the tumor type ($abbreviation)
#' @param useLabels Do you want to use labels within your plot (= TRUE) or only show
#' the different domains with different colors (= FALSE)? Default = TRUE.
#' @return ggplot with all datapoints from the reference cohort,
#' color coded by the domain.
#' @export
#' @import umap
#'
plotCohortDomain <- function(dataUMAPList, useLabels = T) {

  if (require("ggrepel") == F) {
    remotes::install_github("fwallis/ggrepel")
  }

dataUMAP <- dataUMAPList$dataUMAP
  dataLogUMAPlabels <- dataUMAP %>% filter(!(duplicated(Domain)))
  umapCohortDomain <- dataUMAP %>%
    ggplot(aes(x = UMAP1,
               y = UMAP2
    )) +
    theme_classic() +
    labs(color = "Domain") +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =0.7),
          legend.key = element_blank()) +
    geom_point(aes(color = Domain),
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
    umapCohortDomain <- umapCohortDomain +
      ggrepel::geom_label_repel(data = dataLogUMAPlabels,
                                aes(color = Domain,
                                    label = Domain),
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

  return(umapCohortDomain)
}
