#' Create synthetic metadata entries
#'
#' This function is designed to create new metadata entries for
#' tumor (sub)types that only have two samples.
#' Additional entries are needed for analysis of variance (ANOVA),
#' as each group requires at least three entries to calculate a mean and a standard deviation for each gene.
#'
#'
#' @param metaDataRef Original metadata file that contains the tumor (sub)type labels.
#' @param classesWith2 A vector containing the tumor (sub)types that currently
#' only have two samples within the metadata file.
#' @import magrittr
#' @return Metadata file now including extra fake entries, with a minimum of
#' three entries per tumor (sub)type.
#'
createExtraMetaData <- function(metaDataRef,
                                classesWith2
) {

  nestedMetaData <- metaDataRef %>% dplyr::filter(class %in% classesWith2) %>%
    dplyr::group_by(class) %>%
    tidyr::nest(.) %>%
    dplyr::mutate(newLines = purrr::map(data, ~ .x[1,]))

  for (line in base::seq(1:base::nrow(nestedMetaData))) {
    newMetaDataPoint <- base::as.data.frame(nestedMetaData$newLines[[line]])
    newMetaDataPoint$class <- nestedMetaData$class[[line]]
    base::rownames(newMetaDataPoint) <- base::paste0("Synthetic", line)

    if(line == 1) {
      newMetaDataDF <- newMetaDataPoint
    } else {
      newMetaDataDF <- base::rbind(newMetaDataDF, newMetaDataPoint)
    }
  }

  newMetaDataDF <- newMetaDataDF[ , base::colnames(metaDataRef)]


  metaDataExtra <- base::rbind(metaDataRef, newMetaDataDF)

  return(metaDataExtra)
}
