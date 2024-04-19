#' Predict tumor (sub)type for new samples with Minority Classifier
#'
#' @param createdModelsMinority R-object containing the generated RF-models ($modelList), the model for the ribodepletion correction ($riboModelList),
#'  the features that were eventually used for the weighted RF within the different folds ($reducedFeaturesList),
#'   the metadata file associated to the reference cohort ($metaData)
#'  and the metadata for the performed run ($metaDataRun).
#' @param countDataNew Matrix containing the RNA-transcript per million data for the new samples to be classified.
#' Patients are in the columns, different genes in the rows.
#' @param outputDir Directory in which you would like to store the R-object containing the results.
#' @param saveModel Do you want to save the resulting predictions in an R object?
#'
#' @return R-object containing the final classifications ($classifications) for the samples,
#' and the probabilities associated to the different classifications ($probability).
#' @export

newPredictionsMinority <- function(createdModelsMinority, countDataNew,
                                   outputDir,
                                   saveModel = T) {
  # Find the predictions for the test data
  # PREPARE DATA
  countDataNew <- apply(countDataNew,2,function(x) (x/sum(x))*1E6)

  countDataNew <- predictRiboCounts(riboModel = createdModelsMinority$riboModelList$riboModel, data = countDataNew)

  dataLogNew <- log(countDataNew + 1) %>% t() %>% as.data.frame()
  dataLogNew <- dataLogNew[ , createdModelsMinority$reducedFeatures, drop = F]

  # Also create test data and specify which biomaterial IDs are in there
  #testSamples <- rownames(dataLogNew)

  result <- predictTest(modelList = createdModelsMinority$modelList, testData = dataLogNew)

  print("Finished with predicting results")


  classificationList <- convertResultToClassification(result = result,
                                                         metaDataRef = createdModelsMinority$metaDataRef,
                                                         addOriginalCall = F)

  if (saveModel == T) {
    #directory <- outputDir
    filename <- paste0(outputDir, "/minorityClassifierResult.rds")
    if (!dir.exists(outputDir)) {
      dir.create(outputDir) }
  saveRDS(classificationList, file = filename)
  }
  return(classificationList)
}
