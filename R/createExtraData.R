#' Create fake count data entries
#'
#' This function is designed to create new count data entries for
#' tumor (sub)types that only have two samples.
#' Additional entries are needed for analysis of variance (ANOVA),
#' as each group requires at least three entries to calculate a mean and a standard deviation for each gene.
#'  These entries are calculated by taking the mean values for each gene from the two available samples,
#'  giving an 'average' sample as the third input for the ANOVA.
#'  These extra samples are only used during the gene selection procedure,
#'  not during the classification.
#'
#' @param countDataRef Dataframe containing the RNA-transcript per million data. Patient samples are in the rows,
#' different genes in the columns.The data has been cleaned by the ribodepletion correction protocol.
#' @param classesWith2 A vector containing the tumor (sub)types that currently
#' only have two samples within the metadata file.
#'
#' @return Count-data dataframe now including extra fake entries, with a minimum of
#' three entries per tumor (sub)type.
#'
createExtraData <- function(countDataRef,
                            classesWith2) {

  subsetCountData <- countDataRef %>% filter(class %in% classesWith2)

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

  countDataExtra <- rbind(countDataRef, newDataDF)

  colnumber <- which(colnames(countDataExtra) == "class")
  countDataExtra <- countDataExtra[, -(colnumber)]

  return(countDataExtra)
}
