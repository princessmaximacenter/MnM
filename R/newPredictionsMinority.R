#' Predict tumor (sub)type for new samples with Minority Classifier
#'
#' @param createdModelsMinority R-object containing the generated RF-models ($modelList), the model for the ribodepletion correction ($riboModelList),
#'  the features that were eventually used for the weighted RF within the different folds ($reducedFeaturesList),
#'   the metadata file associated to the reference cohort ($metaDataRef)
#'  and the metadata for the performed run ($metaDataRun).
#' @param countDataNew Matrix containing the RNA-transcript per million data for the new samples to be classified.
#' Samples are in the columns, different RNA-transcripts in the rows.
#' @param outputDir Directory in which you would like to store the R-object containing the results.
#' @param saveModel Do you want to save the results in an R object?
#' @param correctRibo Do you want to perform a correction for the ribodepletion protocol on your dataset? Default is TRUE.
#' @return R-object containing the final classifications ($classifications) for the samples,
#' the probabilities associated to the different classifications ($probability),
#' the metadata file associated to the reference cohort ($metaDataRef)
#'  and the metadata for the performed run ($metaDataRun).
#' @export

newPredictionsMinority <- function(createdModelsMinority,
                                   countDataNew,
                                   outputDir = paste0("./", format(as.Date(Sys.Date(), "%Y-%m-%d"), "%Y_%m_%d")),
                                   saveModel = T,
                                   correctRibo = T
                                   ) {
  # Find the predictions for the test data
  # PREPARE DATA
  if (saveModel == T) {
    #directory <- outputDir

    if (!dir.exists(outputDir)) {
      checkDirectory <- base::tryCatch(base::dir.create(outputDir))
      if (checkDirectory == F) {
        stop("The directory you want the classification to be saved in cannot be created due to an error in the directory path. Please check the spelling of your specified outputDir.")
      }
    }
  }

  countDataNew <- base::apply(countDataNew,2,function(x) (x/base::sum(x))*1E6)
  if (correctRibo == T) {
  countDataNew <- predictRiboCounts(riboModel = createdModelsMinority$riboModelList$riboModel, data = countDataNew)
  }
  dataLogNew <- base::log(countDataNew + 1) %>% base::t() %>% base::as.data.frame()
  dataLogNew <- dataLogNew[ , createdModelsMinority$reducedFeatures, drop = F]

  # Also create test data and specify which biomaterial IDs are in there
  #testSamples <- rownames(dataLogNew)

  result <- predictTest(modelList = createdModelsMinority$modelList,
                        testData = dataLogNew)

  base::print("Finished with classifying results")


  classificationList <- convertResultToClassification(result = result,
                                                         metaDataRef = createdModelsMinority$metaDataRef,
                                                         addOriginalCall = F)

  classificationList$metaDataRun <- createdModelsMinority$metaDataRun
  if (saveModel == T) {
    filename <- base::paste0(outputDir, "/minorityClassifierResult.rds")

    base::saveRDS(classificationList, file = filename)
  base::print(base::paste0("Please find the generated R-object with the classification results within ", filename))
  }
  return(classificationList)
}
