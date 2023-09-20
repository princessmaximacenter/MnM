#'  Plot UMAP Domains
#'
#' Function to plot the UMAP coordinates for the different domains.
#' All domains are colored differently to
#' distinguish them from one another.
#'
#' @param dataUMAP Dataframe containing the UMAP coordinates, the tumor type ($subclass),
#' the domain ($Domain) and the desired abbreviation for the tumor type ($abbreviation)
#' @return ggplot with all datapoints from the reference cohort,
#' color coded by the domain.
#'
plotCohortDomain <- function(dataUMAP) {

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
          plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm")) +

    scale_color_manual(values = c("Hemato" = "#880808",
                                  "Solid" =  "#D1944A",
                                  "Neuro" = "#012695"),
                       labels = c("Blood tumors", 'Neurological tumors', "Solid tumors"))

  return(umapCohortDomain)
}
