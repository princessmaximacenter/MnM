#' Plot pie chart from one domain
#'
#' This function is dedicated to plotting the pie chart from one domain.
#'
#' @param data Dataframe containing the fractions, angles and positions of the different labels,
#' with the labels being present on both the classColumn and higherClassColumn level.
#' @param domain Which domain do you want to plot? This name should be present in the
#' domainColumn of the metadata.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param domainColumn Column in the metadata file that contains the domain labels.
#' @param textSizeClass What size should the letters have for the tumor type?
#' @param textSizeSubspec What size should the letters have for the tumor subtype?
#' @param freqSameTumorType How many subtypes belong to one tumor type? In order of the tumor types.
#' @param plotColors Which colors do you want to use to color your plot?
#' @param saveImage Do you want to save your image?
#' @param storeLocation Where do you want to save the image?
#' @param includeNumbers Do you want to show the numbers for the pie chart on the outer edge?
#' @return Piechart of how many samples are present within the dataset belonging to a
#' certain tumor type, and to a certain tumor subtype.
#' The piechart is generated for one domain, as there are many different tumor types and subtypes
#' within each domain.
#' @import grDevices

plotPieChartDomain <- function(data,
                               domain,
                               classColumn,
                               higherClassColumn,
                               domainColumn,
                               textSizeClass,
                               textSizeSubspec,
                               freqSameTumorType,
                               plotColors,
                               saveImage,
                               storeLocation,
                               includeNumbers) {


  for (i in seq(1:length(freqSameTumorType))) {
    extraCols <- base::rep(plotColors[i], freqSameTumorType[i])
    if (i == 1) {
      myColours <- extraCols
    } else {
      myColours <- c(myColours, extraCols)
    }
  }
  filename <- storeLocation
  if (saveImage == T) {
    grDevices::pdf(file = filename ) }
  graphics::par(mar=c(0,0,0,0))
  if (includeNumbers == T) {
  tumorTypeSubtypePie(data %>% dplyr::group_by(!!dplyr::sym(classColumn)) %>%
                        dplyr::select(counts) %>%
                        tibble::deframe(),#sapply(unique(data[,classColumn]), function(x) sum(data$counts[data[,classColumn] == x])),
       init.angle = 0,
       labels = data$counts,
       border = F,
       radius=0.9,
       cex=0.5,
       col=NULL)
  graphics::par(new=T)}

  tumorTypeSubtypePie(data %>% dplyr::group_by(!!dplyr::sym(classColumn)) %>%
                        dplyr::select(counts) %>%
                        tibble::deframe(),
       init.angle = 0,
       labels = base::rep("",base::length(data[,classColumn])),
       radius=0.95,
       col=myColours,
       border = "white")

  graphics::par(new=T)

  tumorTypeSubtypePie(data %>% dplyr::group_by(!!dplyr::sym(classColumn)) %>%
                        dplyr::select(counts) %>%
                        tibble::deframe(), # 2
       init.angle = 0,
       radius=0.38,
       border = F,
       cex=textSizeSubspec,
       col=NULL)
  graphics::par(new=T)

  tumorTypeSubtypePie(sapply(base::unique(data[,higherClassColumn]), function(x) sum(data$counts[data[,higherClassColumn] == x])), # White border around subclass
       init.angle = 0,
       labels = base::rep("",base::length(base::unique(data[,higherClassColumn]))),
       border = "white",radius=0.38,
       col=plotColors)
  graphics::par(new=T)

  tumorTypeSubtypePie(sapply(base::unique(data[,higherClassColumn]),
              function(x) sum(data$counts[data[,higherClassColumn] == x])),
       init.angle = 0,
       radius=0.25,
       border = F,
       cex=textSizeClass,
       col=NULL) # 4
  graphics::par(new=T)

  tumorTypeSubtypePie(sapply(base::unique(data[,domainColumn]),
              function(x) sum(data$counts[data[,domainColumn] == x])),
       init.angle = 0,
       labels = rep("",3),
       border = F,
       radius=0.25,
       col="white")
  graphics::par(new=T)
  graphics::text(0,0,labels = domain,cex=1, col = "grey")
  if (saveImage == T) {
    grDevices::dev.off()
    base::print(base::paste0("Image has been saved under ", filename))
  }
}
