#' Title
#'
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param abbreviationSubtype Dataframe containing the link between the
#' tumor subtype label and tumor subtype abbreviation.
#' @param abbreviationTumorType Dataframe containing the link between the
#' tumor type label and tumor type abbreviation.
#' @param highestClassColumn Potentially not needed.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#'
#' @return Dataframe containing ...
#'
getPieDF <- function(metaDataRef,
                     abbreviationSubtype,
                     abbreviationTumorType,
                     highestClassColumn = "Disease_main_class",
                     higherClassColumn,
                     classColumn
                     ) {
  pieDF <- metaDataRef[!duplicated(metaDataRef[, classColumn]),]
  pieDF$counts <- table(metaDataRef[, classColumn])[pieDF[, classColumn]]
  pieDF$fraction <- pieDF$counts/sum(pieDF$counts)

  #pieDF <- pieDF[order(pieDF$Domain,pieDF[,highestClassColumn],pieDF[,higherClassColumn],pieDF[, classColumn]),]

  pieDF$ymax <- cumsum(pieDF$fraction)

  pieDF$ymin <- c(0,cumsum(pieDF$fraction)[-nrow(pieDF)])


  #pieDF$abbreviation <- NA
  for (i in seq(1:nrow(pieDF))) {
    pieDF[,higherClassColumn][i] <- abbreviationTumorType[abbreviationTumorType[,higherClassColumn] == pieDF[,higherClassColumn][i],"abbreviation"]
    pieDF[,classColumn][i] <- abbreviationSubtype[abbreviationSubtype[,classColumn] == pieDF[,classColumn][i],"abbreviation"]
  }

  return(pieDF)
}
