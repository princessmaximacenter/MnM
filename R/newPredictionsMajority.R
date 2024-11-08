#' Predict tumor subtype for new samples with Majority classifier
#'
#' Function to generate predictions for new input samples from their
#' RNA-seq count data.
#'
#' @param createdModelsMajority  R-object containing the rotations and scalings
#' for each training data subset ($rotationsAndScalingList),
#' the model to correct for the ribodepletion efficacy ($riboModelList),
#' which samples were present in each subset ($samplesTrainDefList),
#' which RNA-transcripts were considered for transformation in the analysis ($nonZeroGenes),
#' the metadata file associated to the reference cohort ($metaDataRef),
#' and metadata for the performed run ($metaDataRun).
#' @param countDataNew Matrix containing the RNA-transcript per million data for the new samples to be classified.
#' Samples are in the columns, different RNA-transcripts in the rows.
#' @param outputDir Directory in which you would like to store the R-object containing the results. Default is today's date.
#' @param saveModel Do you want to save the resulting predictions in an R object? Default is TRUE.
#' @param correctRibo Do you want to perform a correction for the ribodepletion protocol on your dataset? Default is TRUE.
#' @param CPMcorrect Do you want to make sure you have CPM counts? Put false in case when you supply a subset of RNA-transcripts at the start.
#' @param whichKimputation The number of neighbor datapoints that need to be used to calculate missing RNA-transcripts from (in case needed).
#'
#' @return R-object containing the final classifications ($classifications) for the samples,
#' and the probabilities associated to the different classifications ($probability).
#' @export
#'
newPredictionsMajority <- function(createdModelsMajority,
                                   countDataNew,
                                   outputDir = paste0("./", format(as.Date(Sys.Date(), "%Y-%m-%d"), "%Y_%m_%d")),
                                   saveModel = T,
                                   correctRibo = T,
                                   CPMcorrect = T,
                                   whichKimputation = 3
) {
  `%notin%` <<- Negate(`%in%`)
  checkFormatTestData(countDataNew = countDataNew,
                      countDataRef = createdModelsMajority$countDataRef,
                      outputDir = outputDir,
                      saveModel = saveModel)

  # Make sure you have CPM counts
  if (CPMcorrect == T) {
    countDataNew <- base::apply(countDataNew,2,function(x) (x/base::sum(x))*1E6)
  }

  # Here I need to include missing gene imputation

  if (correctRibo == T) {
    base::cat("\nStarting the ribocorrection procedure.\n")
  countDataNew <- predictRiboCounts(riboModel = createdModelsMajority$riboModelList$riboModel,
                                    data = countDataNew,
                                    countDataRef = createdModelsMajority$countDataRef,
                                    whichKimputation = whichKimputation)

  }

  neededGenes <- purrr::map(createdModelsMajority$rotationsAndScalingsList[["scaleFeaturesList"]],
                                          ~.x$varFeatures) %>%
    purrr::reduce( function(x,y) {dplyr::union(x,y)})
  missingGenes <- neededGenes[neededGenes %notin% rownames(countDataNew)]

  if (base::length(missingGenes) > 0) {
    base::cat(base::paste0("There are ",
                           base::length(missingGenes),
                           " genes missing from the dataset.\nImputing their values\n"))
    countDataNew <- calculateMissingGenes(countDataNew = countDataNew,
                                          neededGenes = neededGenes,
                                          countDataRef = createdModelsMajority$riboModelList$counts,
                                          whichK = whichKimputation)
  }


  # Log-transform data
  dataLogNew <- base::log(countDataNew + 1)

  testSamples <- base::colnames(dataLogNew)

  result <- obtainPredictionMajorityClassifier(rotationsAndScalingsList = createdModelsMajority$rotationsAndScalingsList,
                                     dataTest = dataLogNew,
                                     metaDataRef = createdModelsMajority$metaDataRef,
                                     classColumn = createdModelsMajority$metaDataRun$classColumn,
                                     nModels = createdModelsMajority$metaDataRun$nModels,
                                     testSamples = testSamples,
                                     vectorK = createdModelsMajority$vectorK
  )

   classificationList <- convertResultToClassification(result = result,
                                                         metaDataRef = createdModelsMajority$metaDataRef,
                                                         addOriginalCall = F)

   classificationList$metaDataRun <- createdModelsMajority$metaDataRun
if (saveModel == T) {

  filename <- base::paste0(outputDir, "/majorityClassifierResult.rds")
  base::saveRDS(classificationList, file = filename)
  base::cat(base::paste0("Please find the generated R-object with the classification results within ", filename))
}
  return(classificationList)
}
