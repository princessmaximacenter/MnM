#' Title
#'
#' @param metaData
#' @param classesWith2
#'
#' @return
#' @export
#'
#' @examples
createExtraMetaData <- function(metaData,
                                classesWith2
) {

  nestedMetaData <- metaData %>% filter(class %in% classesWith2) %>%
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

  col_order <- colnames(metaData)
  newMetaDataDF <- newMetaDataDF[ , col_order]


  metaDataExtra <- rbind(metaData, newMetaDataDF)

  return(metaDataExtra)
}
