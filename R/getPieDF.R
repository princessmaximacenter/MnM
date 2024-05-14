#' Convert metadata into frequency dataframe
#'
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param abbreviations Dataframe containing the different tumor subtype labels (classColumn) with their appropriate abbreviations ($abbreviationSubtype),
#' and tumor type labels (higherClassColumn) with their appropriate abbreviations (abbreviationTumorType).
#' Please make sure that the higherClassColumn name is used as column name within abbreviationsCombi for the tumor type labels,
#' and classColumn name for the tumor subtype labels.
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
  pieDF <- metaDataRef[!base::duplicated(metaDataRef[, c(classColumn,higherClassColumn)]),]
  pieDF$counts <- NA
  for (tumorType in base::unique(metaDataRef[,higherClassColumn])) {
    pieDF[pieDF[,higherClassColumn] == tumorType,"counts"] <- base::table(metaDataRef[metaDataRef[,higherClassColumn] == tumorType, classColumn])[pieDF[pieDF[,higherClassColumn] == tumorType, classColumn]]

  }
  pieDF$fraction <- pieDF$counts/base::sum(pieDF$counts)

  #pieDF <- pieDF[order(pieDF$Domain,pieDF[,highestClassColumn],pieDF[,higherClassColumn],pieDF[, classColumn]),]

  pieDF$ymax <- base::cumsum(pieDF$fraction)

  pieDF$ymin <- c(0,base::cumsum(pieDF$fraction)[-base::nrow(pieDF)])


  #pieDF$abbreviation <- NA
  for (i in base::seq(1:base::nrow(pieDF))) {
    pieDF[,classColumn][i] <- abbreviations[(abbreviations[,classColumn] == pieDF[,classColumn][i]) & (abbreviations[,higherClassColumn] == pieDF[,higherClassColumn][i]),"abbreviationSubtype"]
    pieDF[,higherClassColumn][i] <- abbreviations[abbreviations[,higherClassColumn] == pieDF[,higherClassColumn][i],"abbreviationTumorType"] %>% base::unique()

  }

  return(pieDF)
}
