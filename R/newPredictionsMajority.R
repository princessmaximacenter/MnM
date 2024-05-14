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
#' @param countDataRef Matrix containing the RNA-transcript per million data for the reference cohort.
#' Samples are in the columns, different RNA-transcripts in the rows.
#' @param countDataNew Matrix containing the RNA-transcript per million data for the new samples to be classified.
#' Samples are in the columns, different RNA-transcripts in the rows.
#' @param outputDir Directory in which you would like to store the R-object containing the results. Default is today's date.
#' @param saveModel Do you want to save the resulting predictions in an R object? Default is TRUE.
#'
#' @return R-object containing the final classifications ($classifications) for the samples,
#' and the probabilities associated to the different classifications ($probability).
#' @export
#'
newPredictionsMajority <- function(createdModelsMajority,
                                   countDataRef,
                                   countDataNew,
                                   outputDir = paste0("./", format(as.Date(Sys.Date(), "%Y-%m-%d"), "%Y_%m_%d")),
                                   saveModel = T
) {

  if (saveModel == T) {
    #directory <- outputDir

    if (!dir.exists(outputDir)) {
      checkDirectory <- base::tryCatch(base::dir.create(outputDir))
      if (checkDirectory == F) {
        stop("The directory you want the classification to be saved in cannot be created due to an error in the directory path. Please check the spelling of your specified outputDir.")
      }
    }
  }

  # Make sure you have CPM counts
  countDataNew <- apply(countDataNew,2,function(x) (x/sum(x))*1E6)
  countDataRef <- apply(countDataRef,2,function(x) (x/sum(x))*1E6)


  countDataNew <- predictRiboCounts(riboModel = createdModelsMajority$riboModelList$riboModel,
                                    data = countDataNew)
  countDataRef<- predictRiboCounts(riboModel = createdModelsMajority$riboModelList$riboModel,
                                   data = countDataRef)


  # Log-transform data
  dataLogRef <- log(countDataRef +1)
  dataLogNew <- log(countDataNew + 1)

  dataLogNonZero <- dataLogRef[createdModelsMajority$nonZeroGenes,]
  dataLogNewNonZero <- dataLogNew[createdModelsMajority$nonZeroGenes, , drop = F]
  testSamples <- colnames(dataLogNewNonZero)

  result <- obtainPredictionMajorityClassifier(rotationsAndScalingsList = createdModelsMajority$rotationsAndScalingsList,
                                     dataTrain = dataLogNonZero,
                                     dataTest = dataLogNewNonZero,
                                     metaDataRef = createdModelsMajority$metaDataRef,
                                     classColumn = createdModelsMajority$metaDataRun$classColumn,
                                     nModels = createdModelsMajority$metaDataRun$nModels,
                                     testSamples = testSamples,
                                     maxNeighbors = createdModelsMajority$metaDataRun$maxNeighbors
  )

   classificationList <- convertResultToClassification(result = result,
                                                         metaDataRef = createdModelsMajority$metaDataRef,
                                                         addOriginalCall = F)

   classificationList$metaDataRun <- createdModelsMajority$metaDataRun
if (saveModel == T) {

  if (!dir.exists(outputDir)) {
    dir.create(outputDir) }

  filename <- paste0(outputDir, "/majorityClassifierResult.rds")
  saveRDS(classificationList, file = filename)
  print(paste0("Please find the generated R-object with the classification results within ", filename))
}
  return(classificationList)
}
