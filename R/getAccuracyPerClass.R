#' Obtain accuracy per tumor (sub)type
#'
#' @param predictions Dataframe containing the diagnosis label ($originalCall) and the top-ranked classification label ($predict) from a classifier.
#'
#' @return Dataframe containing the accuracy for the different tumor (sub)types within the dataset.
#' @export
getAccuracyPerClass <- function(predictions) {

  classes <- predictions[, "originalCall"] %>% base::unique()


  for (i in base::seq(1:base::length(classes))) {
    predictionsFiltered <- predictions %>% dplyr::filter(originalCall == classes[i])
    correct <- predictionsFiltered %>% dplyr::filter(originalCall == predict)

    accuracyDF <- base::data.frame(
      accuracy = base::nrow(correct) / base::nrow(predictionsFiltered))

    base::rownames(accuracyDF) <- classes[i]
    if (i == 1) {
      accuracyDFTotal <- accuracyDF
    } else {
      accuracyDFTotal <- base::rbind(accuracyDFTotal,
                               accuracyDF)
    }
  }
  return(accuracyDFTotal)

}
