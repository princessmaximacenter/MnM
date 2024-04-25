#' Obtain accuracy per tumor (sub)type
#'
#' @param predictions Dataframe containing the diagnosis label ($originalCall) and the top-ranked classification label ($predict) from a classifier.
#'
#' @return Dataframe containing the accuracy for the different tumor (sub)types within the dataset.
#' @export
getAccuracyPerClass <- function(predictions) {

  classes <- predictions[, "originalCall"] %>% unique()


  for (i in seq(1:length(classes))) {
    predictionsFiltered <- predictions %>% dplyr::filter(originalCall == classes[i])
    correct <- predictionsFiltered %>% dplyr::filter(originalCall == predict)

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
