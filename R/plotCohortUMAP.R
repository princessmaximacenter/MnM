#' Plot UMAP coordinates colored by tumor domain, tumor types or tumor subtypes
#'
#' Function to plot the UMAP coordinates for the different tumor (sub)type samples
#' belonging to one domain or tumor type.
#'
#' @param dataUMAPList List containing the UMAP coordinate data
#' and the desired abbreviation for the tumor type ($abbreviation), resulting from the function 'createUMAPcohort'
#' @param domain Which domain do you want to plot the tumor types for? If not specified, visualization will be performed on the domain level by default.
#' @param subtype In case domain is specified: Do you want to visualize on the tumor subtype level (TRUE) or on the tumor type level (FALSE)?
#' @param tumorType If domain is specified and subtype = TRUE, there's a possibility to plot the tumor subtypes for only one tumor type.
#' Please specify which tumor type you would like to visualize within this parameter.
#' @param plotColors Which colors do you want to use for the tumor types? If not specified, default ggplot2 colors will be used.
#' @param useLabels Do you want to supply labels within the plot?

#' If that would be desired, fill out the tumor type that you would like to visualize at 'tumorType' (for example: "B-ALL").

#' @return ggplot with only the datapoints of the selected domain,
#' color coded by the tumor type (subtype = F) or tumor subtype (subtype = T).
#' @export
#' @import remotes
#'
plotCohortUMAP <- function(dataUMAPList,
                                domain = NA,
                                subtype = F,
                                tumorType = NA,
                                plotColors = NA,
                                useLabels = T
                                ) {

  if (require("ggrepel") == F) {
    remotes::install_github("fwallis/ggrepel")
  }

  dataUMAP <- dataUMAPList$dataUMAP

  if (!base::is.na(domain)) {
    umapDomain <- dataUMAP %>% dplyr::filter(Domain == domain)
    abbreviations <- dataUMAPList$abbreviations %>% dplyr::filter(abbreviationSubtype %in% base::unique(umapDomain$abbreviationSubtype))

    if (subtype == F) {
      umapDomain$abbreviation <- umapDomain$abbreviationTumorType
      domainTumorTypes <- abbreviations %>%
        dplyr::filter(Domain == domain) %>%
        dplyr::select(abbreviationTumorType) %>%
        base::unique() %>%
        tibble::deframe()
    } else {
      umapDomain$abbreviation <- umapDomain$abbreviationSubtype
      if (!base::is.na(tumorType)[1]) {
        umapDomain <- umapDomain %>% dplyr::filter(subclass == tumorType)
      }
      domainTumorTypes <- abbreviations %>%
        dplyr::filter(Domain == domain) %>%
        dplyr::select(abbreviationSubtype) %>%
        base::unique() %>%
        tibble::deframe()

    }
  } else {
    umapDomain <- dataUMAP
    umapDomain$abbreviation <- umapDomain$Domain
  }

  dataLogUMAPlabels <- umapDomain %>% dplyr::filter(!(duplicated(abbreviation)))

  dataLogUMAPlabels %<>% dplyr::arrange(subclass)
  dataLogUMAPlabels$abbreviation <- base::factor(dataLogUMAPlabels$abbreviation,
                                                 levels = base::unique(dataLogUMAPlabels$abbreviation))


  umapCohortTumorType <- umapDomain %>%
    ggplot2::ggplot(
      ggplot2::aes(x = UMAP1,
               y = UMAP2
    )) +
    ggplot2::theme_classic() +
    ggplot2::labs(color = "Tumor Type") +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(colour = "black", fill=NA, linewidth =0.7),
          legend.key = ggplot2::element_blank(),
          legend.background = ggplot2::element_blank(),
          axis.title = ggplot2::element_text(size = 14),
          axis.title.x = ggplot2::element_text(vjust = -1.8),
          axis.title.y = ggplot2::element_text(vjust = 2),
          legend.position = "none",
          plot.margin = ggplot2::unit(c(0.8,0.8,0.8,0.8), "cm")) +
    ggplot2::geom_point(ggplot2::aes(color = abbreviation),
               shape = 19,
               key_glyph = ggplot2::draw_key_point)

  if (useLabels == T) {
    umapCohortTumorType <- umapCohortTumorType +
      ggrepel::geom_label_repel(data = dataLogUMAPlabels,
                                ggplot2::aes(color = abbreviation,
                           label = abbreviation),
                       max.overlaps = 40,
                       size=4,
                       seed = 1,
                       label.size = 1,
                       show.legend = F,
                       fill = NA,
                       segment.alpha = 0
      )
  }

  if (!base::is.na(plotColors)[1]) {
    umapCohortTumorType <- umapCohortTumorType +
      ggplot2::scale_color_manual(values = plotColors,
                         breaks = domainTumorTypes)
  }

  return(umapCohortTumorType)


}
