#' Create fake metadata entries
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

  nestedMetaData <- metaDataRef %>% filter(class %in% classesWith2) %>%
    group_by(class) %>%
    nest(.) %>%
    mutate(newLines = map(data, ~ .x[1,]))

  for (line in seq(1:nrow(nestedMetaData))) {
    newMetaDataPoint <- as.data.frame(nestedMetaData$newLines[[line]])
    newMetaDataPoint$class <- nestedMetaData$class[[line]]
    rownames(newMetaDataPoint) <- paste0("Fake", line)

    if(line == 1) {
      newMetaDataDF <- newMetaDataPoint
    } else {
      newMetaDataDF <- rbind(newMetaDataDF, newMetaDataPoint)
    }
  }

  col_order <- colnames(metaDataRef)
  newMetaDataDF <- newMetaDataDF[ , col_order]


  metaDataExtra <- rbind(metaDataRef, newMetaDataDF)

  return(metaDataExtra)
}
