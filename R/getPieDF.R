#' Convert metadata into frequency dataframe
#'
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param abbreviations Dataframe containing the different tumor subtype labels with their appropriate abbreviations ($abbreviationSubtype),
#' and tumor type labels with their appropriate abbreviations (abbreviationTumorType).
#' Please make sure that the higherClassColumn name is used as column name within abbreviationsCombi for the tumor type labels,
#' and classColumn name for the tumor subtype labels.
#' tumor type label and tumor type abbreviation.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#'
#' @return Dataframe containing the different labels within our dataset with their frequencies.
#'
getPieDF <- function(metaDataRef,
                     abbreviations,
                     higherClassColumn,
                     classColumn
                     ) {
  pieDF <- metaDataRef[!duplicated(metaDataRef[, c(classColumn,higherClassColumn)]),]
  #pieDF$counts <- table(metaDataRef[, classColumn])[pieDF[, classColumn]]
  pieDF$counts <- NA
  for (sth in unique(metaDataRef[,higherClassColumn])) {
    pieDF[pieDF[,higherClassColumn] == sth,"counts"] <- table(metaDataRef[metaDataRef[,higherClassColumn] == sth, classColumn])[pieDF[pieDF[,higherClassColumn] == sth, classColumn]]

  }
  pieDF$fraction <- pieDF$counts/sum(pieDF$counts)

  #pieDF <- pieDF[order(pieDF$Domain,pieDF[,highestClassColumn],pieDF[,higherClassColumn],pieDF[, classColumn]),]

  pieDF$ymax <- cumsum(pieDF$fraction)

  pieDF$ymin <- c(0,cumsum(pieDF$fraction)[-nrow(pieDF)])


  #pieDF$abbreviation <- NA
  for (i in seq(1:nrow(pieDF))) {
    pieDF[,classColumn][i] <- abbreviations[(abbreviations[,classColumn] == pieDF[,classColumn][i]) & (abbreviations[,higherClassColumn] == pieDF[,higherClassColumn][i]),"abbreviationSubtype"]
    pieDF[,higherClassColumn][i] <- abbreviations[abbreviations[,higherClassColumn] == pieDF[,higherClassColumn][i],"abbreviationTumorType"] %>% unique()

  }

  return(pieDF)
}
