#' Plot pie-chart per domain as cohort overview
#'
#' Function to automatically generate a pie-chart for the samples that belong to
#' one of the tumor domains within the dataset. Required are dedicated dataframes
#' to specify the names that you would like to use for the different tumor types and subtypes,
#' as usually their full names are too lengthy to fit on the plot.
#'
#' Also, the colors you want to use within the plot could be specified,
#' so that the sheer amount of tumor types does not cause a decreased visibility.
#'
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param whichDomain Which domain do you want to plot? This name should be present in the
#' domainColumn of the metadata.
#' @param storeLocation Where do you want to store your generated image? Default is within working directory.
#' @param saveImage Do you want to save your image? Default is FALSE.
#' @param textSizeClass What size should the letters have for the tumor type? Default is 0.5.
#' @param textSizeSubspec What size should the letters have for the tumor subtype? Default is 0.45.
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Name of column in the metadata file that contains the tumor type labels.
#' @param domainColumn Name of column in the metadata file that contains the domain labels.
#' @param abbreviations Dataframe containing the different tumor subtype labels with their appropriate abbreviations ($abbreviationSubtype),
#' and tumor type labels with their appropriate abbreviations (abbreviationTumorType).
#' Please make sure that the higherClassColumn name is used as column name within abbreviations for the tumor type labels,
#' and classColumn name for the tumor subtype labels.
#' @param plotColors Which colors do you want to use to color your plot? If not specified,
#' default colors from RColorBrewer palette spectral will be used.
#' @param includeNumbers  Do you want to show the numbers for the pie chart on the outer edge?
#' Default is includeNumbers = TRUE.
#' @param lowerLimitFraction Tumor subtype fraction that is required for plotting the
#' name on the pie-chart. In case the fraction is lower than the limit, the name will be removed.
#'
#' @return Pie-chart of how many samples are present within the dataset belonging to a
#' certain tumor type, and to a certain tumor subtype.
#' The pie-chart is generated for one domain, as there are many different tumor types and subtypes
#' within each domain.
#' @export
#' @import RColorBrewer common
#'
createPieChartImage <- function(metaDataRef,
         classColumn,
         higherClassColumn,
         domainColumn = NA,
         abbreviations = NA,
         whichDomain = NA,
         storeLocation = ".",
         saveImage = F,
         textSizeClass = 0.5,
         textSizeSubspec = 0.45,
         plotColors = NA,
         includeNumbers = T,
         lowerLimitFraction = 0
         ) {


  if (is.na(abbreviations)[1] & !is.na(domainColumn)[1]) {

      abbreviations <- metaDataRef[, c(domainColumn, classColumn, higherClassColumn)] %>% base::unique()
      abbreviations$abbreviationSubtype <- abbreviations[,classColumn]
      abbreviations$abbreviationTumorType <- abbreviations[,higherClassColumn]

      print("You have not supplied any abbreviations (abbreviations). If you would like to use abbreviations,
          please generate a dataframe with the following columns and abbreviations within abbreviationSubtype and abbreviationTumorType:")
      print(abbreviations[1:4,])
  } else if (is.na(abbreviations)[1] ) {
    abbreviations <- metaDataRef[, c(classColumn, higherClassColumn)] %>% base::unique()
    abbreviations$domainColumn <- " "
    abbreviations$abbreviationSubtype <- abbreviations[,classColumn]
    abbreviations$abbreviationTumorType <- abbreviations[,higherClassColumn]

    print("You have not supplied any abbreviations (abbreviations). If you would like to use abbreviations,
          please generate a dataframe with the following columns and abbreviations within abbreviationSubtype and abbreviationTumorType:")
    print(abbreviations[1:4,])
  }


  if (is.na(domainColumn)[1]) {
    domainColumn <- "domainColumn"
    whichDomain <- " "
    metaDataRef$domainColumn <- whichDomain
    }

  abbreviationsToRemove <- pieDF %>% dplyr::filter(fraction < lowerLimitFraction) %>% dplyr::pull(site_group)
  for (i in seq_along(abbreviationsToRemove)) {
    abbreviations[abbreviations$site_group == abbreviationsToRemove[i],"abbreviationSubtype"] <- common::spaces(i)
  }

  pieDF <- getPieDF(metaDataRef = metaDataRef,
                    abbreviations = abbreviations,
                    classColumn = classColumn,
                    higherClassColumn = higherClassColumn
  )

  pieDF <- pieDF %>% dplyr::arrange(!!sym(higherClassColumn), !!sym(classColumn))

  abbreviations %<>% arrange(!!sym(higherClassColumn), !!sym(classColumn))

 # if(!is.na(domainColumn)[1]) {
 #   abbreviationCombi %<>% dplyr::filter(!!dplyr::sym(domainColumn) == whichDomain)
#  }

namesDomain <- abbreviations %>% dplyr::filter(!!dplyr::sym(domainColumn) == whichDomain,
                                               !!dplyr::sym(higherClassColumn) %in% base::unique(metaDataRef[,higherClassColumn])) %>%
  dplyr::select(abbreviationTumorType) %>% base::unique() %>%
  tibble::deframe()

namesSubtype <- abbreviations %>% dplyr::filter(!!dplyr::sym(domainColumn) == whichDomain,
                                                !!dplyr::sym(classColumn) %in% base::unique(metaDataRef[,classColumn])) %>%
  dplyr::select(abbreviationSubtype) %>% base::unique() %>%
  tibble::deframe()
pieChartDF <- pieDF[pieDF[,domainColumn] == whichDomain,]
#pieChartDF <- pieChartDF %>% dplyr::slice(order(factor(!!dplyr::sym(higherClassColumn), levels = namesDomain)))
pieChartDF <- pieChartDF %>% dplyr::slice(order(factor(!!dplyr::sym(classColumn),
                                                       levels = namesSubtype)))

freqSameTumorType <- pieChartDF[,higherClassColumn] %>% base::table()
freqSameTumorType <- freqSameTumorType[base::match(namesDomain, base::names(freqSameTumorType))]

if (is.na(plotColors)[1]) {
  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(10, "Spectral"))

  plotColors <- getPalette(base::length(freqSameTumorType))
}

plotPieChartDomain(data = pieChartDF,
                     domain = whichDomain,
                   classColumn = classColumn,
                   higherClassColumn = higherClassColumn,
                   domainColumn = domainColumn,
                     textSizeClass = textSizeClass,
                     textSizeSubspec = textSizeSubspec,
                   freqSameTumorType = freqSameTumorType,
                   plotColors = plotColors,
                     saveImage = saveImage,
                    storeLocation = storeLocation,
                   includeNumbers = includeNumbers

)
}
