#' Create synthetic count data entries
#'
#' This function is designed to create new count data entries for
#' tumor subtypes that only have two samples.
#' Additional entries are needed for analysis of variance (ANOVA),
#' as each group requires at least three entries for each RNA-transcript.
#'  These synthetic count data entries are calculated by taking the mean values for each RNA-transcript from the two available samples,
#'  giving an 'average' sample as the third input for the ANOVA.
#'  These extra samples are only used during the feature selection procedure,
#'  not during the classification.
#'
#' @param countDataRef Dataframe containing the RNA-transcript per million data. Patient samples are in the rows,
#' different RNA-transcripts in the columns.The data has been cleaned by the ribodepletion correction protocol.
#' @param classesWith2 A vector containing the tumor subtypes that currently
#' only have two samples within the metadata file.
#'
#' @return Count dataframe now including extra synthetic data entries, guaranteeing a minimum of
#' three entries per tumor subtype.
#'
createExtraData <- function(countDataRef,
                            classesWith2) {

  subsetCountData <- countDataRef %>% dplyr::filter(class %in% classesWith2)

  createData <- subsetCountData %>%
    dplyr::group_by(class) %>%
    tidyr::nest(.) %>%
    dplyr::mutate(meanDataPoint = purrr::map(data, ~ base::apply(.x, 2, base::mean)))

  for (line in base::seq(1:base::nrow(createData))) {

    newDataPoint <- base::data.frame(createData$meanDataPoint[[line]])
    base::colnames(newDataPoint) <- createData$class[line]
    if(line == 1) {
      newDataDF <- newDataPoint
    } else {
      newDataDF <- base::cbind(newDataDF, newDataPoint)
    }

  }
  newDataDF <- base::t(newDataDF) %>% base::as.data.frame()
  newDataDF$class <- base::rownames(newDataDF)
  base::rownames(newDataDF) <- base::paste0("Synthetic", 1:base::length(newDataDF$class))

  countDataExtra <- base::rbind(countDataRef, newDataDF)

  colnumber <- base::which(base::colnames(countDataExtra) == "class")
  countDataExtra <- countDataExtra[, -(colnumber)]

  return(countDataExtra)
}
