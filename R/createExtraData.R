#' Title
#'
#' @param countData
#' @param classesWith2
#'
#' @return
#' @export
#'
createExtraData <- function(countData,
                            classesWith2) {

  subsetCountData <- countData %>% filter(class %in% classesWith2)

  createData <- subsetCountData %>%
    group_by(class) %>% nest %>%
    mutate(meanDataPoint = map(data, ~ apply(.x, 2, mean)))

  for (line in seq(1:nrow(createData))) {

    newDataPoint <- data.frame(createData$meanDataPoint[[line]])
    colnames(newDataPoint) <- createData$class[line]
    if(line == 1) {
      newDataDF <- newDataPoint
    } else {
      newDataDF <- cbind(newDataDF, newDataPoint)
    }

  }
  newDataDF <- t(newDataDF) %>% as.data.frame()
  newDataDF$class <- rownames(newDataDF)
  rownames(newDataDF) <- paste0("Fake", 1:length(newDataDF$class))

  countDataExtra <- rbind(countData, newDataDF)

  colnumber <- which(colnames(countDataExtra) == "class")
  countDataExtra <- countDataExtra[, -(colnumber)]

  return(countDataExtra)
}
