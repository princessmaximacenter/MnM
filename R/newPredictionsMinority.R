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
#' @param CPMcorrect Do you want to make sure you have CPM counts? Put false in case when you supply a subset of RNA-transcripts at the start.
#' @return R-object containing the final classifications ($classifications) for the samples,
#' the probabilities associated to the different classifications ($probability),
#' the metadata file associated to the reference cohort ($metaDataRef)
#'  and the metadata for the performed run ($metaDataRun).
#' @export

newPredictionsMinority <- function(createdModelsMinority,
                                   countDataNew,
                                   outputDir = paste0("./", format(as.Date(Sys.Date(), "%Y-%m-%d"), "%Y_%m_%d")),
                                   saveModel = T,
                                   correctRibo = T,
                                   CPMcorrect = T,
                                   whichKimputation = 3
                                   ) {
  # Find the predictions for the test data
  # PREPARE DATA
  `%notin%` <<- Negate(`%in%`)
  checkFormatTestData(countDataNew = countDataNew,
                     countDataRef = createdModelsMinority$countDataRef,
                     outputDir = outputDir,
                     saveModel = saveModel)


  # Make sure you have CPM counts
  if (CPMcorrect == T) {
    countDataNew <- base::apply(countDataNew,2,function(x) (x/base::sum(x))*1E6)
  }



  # Ribo correction if needed
  if (correctRibo == T) {
    cat("\nStarting the ribocorrection procedure.\n")
    countDataNew <- predictRiboCounts(riboModel = createdModelsMinority$riboModelList$riboModel,
                                      data = countDataNew,
                                      countDataRef = createdModelsMinority$countDataRef,
                                      whichKimputation = whichKimputation)

  }

  # Missing gene imputation
  neededGenes <- createdModelsMinority$reducedFeatures
  missingGenes <- neededGenes[neededGenes %notin% rownames(countDataNew)]

  if (base::length(missingGenes) > 0) {
    cat(paste0("There are ", length(missingGenes), " genes missing from the dataset.\nImputing their values.\n"))
    countDataRef <- createdModelsMinority$riboModelList$counts

    countDataNew <- calculateMissingGenes(countDataNew = countDataNew,
                                          neededGenes = neededGenes,
                                          countDataRef = countDataRef,
                                          whichK = whichKimputation)

  }

  # Log transformation
  dataLogNew <- base::log(countDataNew + 1) %>% base::t() %>% base::as.data.frame()
  dataLogNew <- dataLogNew[ , createdModelsMinority$reducedFeatures, drop = F]

  # Get classifications from the individual models for new dataset
  result <- predictTest(modelList = createdModelsMinority$modelList,
                        testData = dataLogNew)

  base::cat("\nFinished with classifying results.")

  # Convert individual model classifications into majority voted classification
  classificationList <- convertResultToClassification(result = result,
                                                         metaDataRef = createdModelsMinority$metaDataRef,
                                                         addOriginalCall = F)

  # Store settings used during the classification procedure
  classificationList$metaDataRun <- createdModelsMinority$metaDataRun
  if (saveModel == T) {
    filename <- base::paste0(outputDir, "/minorityClassifierResult.rds")

    base::saveRDS(classificationList, file = filename)
  base::cat(base::paste0("\nPlease find the generated R-object with the classification results within ", filename))
  }
  return(classificationList)
}
