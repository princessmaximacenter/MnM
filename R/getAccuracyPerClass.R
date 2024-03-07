getAccuracyPerClass <- function(predictions,
                                metaDataRef,
                                classColumn) {

  classes <- metaDataRef[, classColumn] %>% unique()


  for (i in seq(1:length(classes))) {
    predictionsFiltered <- predictions %>% filter(originalCall == classes[i])
    correct <- predictionsFiltered %>% filter(originalCall == predict)

    accuracyDF <- data.frame(
      accuracy = nrow(correct) / nrow(predictionsFiltered))

    rownames(accuracyDF) <- classes[i]
    if (i == 1) {
      accuracyDFTotal <- accuracyDF
    } else {
      accuracyDFTotal <- rbind(accuracyDFTotal,
                               accuracyDF)
    }
  }
  return(accuracyDFTotal)

}
