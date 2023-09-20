#' Plot piechart per domain as cohort overview
#'
#' Function to automatically generate a pie-chart for the samples that belong to
#' one of the tumor domains within the dataset. Required are dedicated dataframes
#' to specify the names that you would like to use for the different tumor (sub)types,
#' as usually their full names are too lengthy to fit on the plot.
#'
#' Also, the colors you want to use within the plot need to be specified,
#' so that the sheer amount of tumor types does not cause a decreased visibility.
#'
#' If you're not sure on which colors to use, check the loadColors function for
#' potential combinations of colors.
#'
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param whichDomain Which domain do you want to plot? This name should be present in the
#' Domain column of the domainColumn of the metadata.
#' @param storeLocation Where do you want to store your generated image?
#' @param saveImage Do you want to save your image? Boolean (T/F) input.
#' @param textSizeClass What size should the letters have for the tumor type?
#' @param textSizeSubspec What size should the letters have for the tumor subtype?
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param domainColumn Column in the metadata file that contains the domain labels.
#' @param abbreviationSubtype Dataframe containing the links between the tumor subtype,
#' the abbreviation required in the plot, the tumor type and the domain.
#' @param abbreviationTumorType Dataframe containing the links between the tumor type,
#' the abbreviation required in the plot, and the domain.
#' @param domainCol Which colors do you want to use to color your plot?
#'
#' @return Piechart of how many samples are present within the dataset belonging to a
#' certain tumor type, and to a certain tumor subtype.
#' The piechart is generated for one domain, as there are many different tumor types and subtypes
#' within each domain.
#' @export
#' @import tibble
#'
createPieChartImage <- function(metaDataRef,
         classColumn,
         higherClassColumn,
         domainColumn,
         abbreviationSubtype,
         abbreviationTumorType,
         whichDomain,
         storeLocation = ".",
         saveImage = F,
         textSizeClass = 0.5,
         textSizeSubspec = 0.5,
         domainCol) {

  pieDF <- getPieDF(metaDataRef = metaDataRef,
                    abbreviationSubtype = abbreviationSubtype,
                    abbreviationTumorType = abbreviationTumorType,
                    classColumn = classColumn,
                    higherClassColumn = higherClassColumn
  )


namesDomain <- abbreviationTumorType %>% filter(!!sym(domainColumn) == whichDomain) %>%
  select(abbreviation) %>%
  deframe()

pieChartDF <- pieDF[pieDF[,domainColumn] == whichDomain,]
pieChartDF <- pieChartDF %>% slice(order(factor(!!sym(higherClassColumn), levels = namesDomain)))

freqSameTumorType <- pieChartDF[,higherClassColumn] %>% table
freqSameTumorType <- freqSameTumorType[match(namesDomain, names(freqSameTumorType))]

plotPieChartDomain(data = pieChartDF,
                     domain = whichDomain,
                   classColumn = classColumn,
                   higherClassColumn = higherClassColumn,
                   domainColumn = domainColumn,
                     textSizeClass = textSizeClass,
                     textSizeSubspec = textSizeSubspec,
                     freq = freqSameTumorType,
                   domainCol = domainCol,
                     saveImage = saveImage,
                    storeLocation = storeLocation
)
}
