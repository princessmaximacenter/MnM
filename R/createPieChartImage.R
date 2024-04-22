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
#' @param abbreviations Dataframe containing the different tumor subtype labels with their appropriate abbreviations ($abbreviationSubtype),
#' and tumor type labels with their appropriate abbreviations (abbreviationTumorType).
#' Please make sure that the higherClassColumn name is used as column name within abbreviations for the tumor type labels,
#' and classColumn name for the tumor subtype labels.
#' @param plotColors Which colors do you want to use to color your plot?
#' @param includeNumbers  Do you want to show the numbers for the pie chart on the outer edge?
#' Default is includeNumbers = TRUE.
#'
#' @return Piechart of how many samples are present within the dataset belonging to a
#' certain tumor type, and to a certain tumor subtype.
#' The piechart is generated for one domain, as there are many different tumor types and subtypes
#' within each domain.
#' @export
#' @import magrittr
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
         textSizeSubspec = 0.5,
         plotColors,
         includeNumbers = T) {


  if (is.na(abbreviations)[1] & !is.na(domainColumn)[1]) {

      abbreviations <- metaDataRef[, c(domainColumn, classColumn, higherClassColumn)] %>% unique()
      abbreviations$abbreviationSubtype <- abbreviations[,classColumn]
      abbreviations$abbreviationTumorType <- abbreviations[,higherClassColumn]

      print("You have not supplied any abbreviations (abbreviations). If you would like to use abbreviations,
          please generate a dataframe with the following columns and abbreviations within abbreviationSubtype and abbreviationTumorType:")
      print(abbreviations[1:4,])
  } else if (is.na(abbreviations)[1] ) {
    abbreviations <- metaDataRef[, c(classColumn, higherClassColumn)] %>% unique()
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



  pieDF <- getPieDF(metaDataRef = metaDataRef,
                    abbreviations = abbreviations,
                    classColumn = classColumn,
                    higherClassColumn = higherClassColumn
  )


 # if(!is.na(domainColumn)[1]) {
 #   abbreviationCombi %<>% dplyr::filter(!!dplyr::sym(domainColumn) == whichDomain)
#  }

namesDomain <- abbreviations %>% dplyr::filter(!!dplyr::sym(domainColumn) == whichDomain) %>%
  dplyr::select(abbreviationTumorType) %>% unique() %>%
  tibble::deframe()

pieChartDF <- pieDF[pieDF[,domainColumn] == whichDomain,]
pieChartDF <- pieChartDF %>% dplyr::slice(order(factor(!!dplyr::sym(higherClassColumn), levels = namesDomain)))

freqSameTumorType <- pieChartDF[,higherClassColumn] %>% table
freqSameTumorType <- freqSameTumorType[match(namesDomain, names(freqSameTumorType))]

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
                    storeLocation = storeLocation

)
}
