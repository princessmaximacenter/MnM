#' Obtain accuracy per tumor (sub)type
#'
#' @param predictions Dataframe containing the diagnosis label ($originalCall) and the top-ranked classification label ($predict) from a classifier.
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#'
#' @return Dataframe containing the accuracy for the different tumor (sub)types within the dataset.
#' @export
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
