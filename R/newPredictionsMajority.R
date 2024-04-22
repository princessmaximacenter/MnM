#' Predict tumor (sub)type for new samples with Majority Classifier
#'
#' Function to generate predictions for new input samples from their
#' RNA-seq count data.
#'
#' @param createdModelsMajority  R-object containing the rotations and scalings
#' for each reference cohort subset($rotationsAndScalingList),
#' the model to correct for the ribodepletion efficacy ($riboModelList),
#' which samples were present in each subset ($samplesTrainDefList),
#' which genes were considered for transformation in the analysis ($nonZeroGenes),
#' the metadata file associated to the reference cohort ($metaData),
#' and metadata for the performed run ($metaDataRun).
#' @param countDataRef Matrix containing the RNA-transcript per million data for the reference cohort.
#' Patients are in the columns,
#' different genes in the rows.
#' @param countDataNew Matrix containing the RNA-transcript per million data for the new samples to be classified.
#' Patients are in the columns, different genes in the rows.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param outputDir Directory in which you would like to store the R-object containing the results.
#' @param saveModel Do you want to save the resulting predictions in an R object?
#'
#' @return R-object containing the final classifications ($classifications) for the samples,
#' and the probabilities associated to the different classifications ($probability).
#' @export
#'
newPredictionsMajority <- function(createdModelsMajority,
                                   countDataRef,
                                   countDataNew,
                                   classColumn,
                                   outputDir = "./",
                                   saveModel = T
) {
  # Make sure you have CPM counts
  countDataNew <- apply(countDataNew,2,function(x) (x/sum(x))*1E6)
  countDataRef <- apply(countDataRef,2,function(x) (x/sum(x))*1E6)


  countDataNew <- predictRiboCounts(riboModel = createdModelsMajority$riboModelList$riboModel, data = countDataNew)
  countDataRef<- predictRiboCounts(riboModel = createdModelsMajority$riboModelList$riboModel, data = countDataRef)


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
                                     samplesTrainDefList = createdModelsMajority$samplesTrainDefList,
                                     classColumn = classColumn,
                                     nModels = createdModelsMajority$metaDataRun$nModels,
                                     testSamples = testSamples,
                                     maxNeighbours = createdModelsMajority$metaDataRun$maxNeighbours
  )

   classificationList <- convertResultToClassification(result = result,
                                                         metaDataRef = createdModelsMajority$metaDataRef,
                                                         addOriginalCall = F)

   classificationList$metaDataRun <- createdModelsMajority$metaDataRun
if (saveModel == T) {
  #directory <- paste0(outputDir, format(as.Date(Sys.Date(), "%Y-%m-%d"), "%m_%d_%Y"))
  if (!dir.exists(outputDir)) {
    dir.create(outputDir) }

  filename <- paste0(outputDir, "/majorityClassifierResult.rds")
  saveRDS(classificationList, file = filename)
}
  return(classificationList)
}
