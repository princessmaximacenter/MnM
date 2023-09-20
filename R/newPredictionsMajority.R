#' Predict tumor (sub)type for new samples
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
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param nModels How many models should be created for the majority voting system?
#' @param nComps How many principal components will be selected after PCA?
#' @param maxNeighbours What is the maximum number of neigbours to be used for the weighted _k_-nearest neighbor algorithm?
#' @param outputDir Directory in which you would like to store the R-object containing the results.
#'
#' @return R-object containing the final classifications ($classifications) for the samples,
#' and the probabilities associated to the different classifications ($probability).
#' @export
#'
newPredictionsMajority <- function(createdModelsMajority = createdModelsMajority,
                                   countDataRef,
                                   countDataNew,
                                   metaDataRef,
                                   classColumn,
                                   nModels,
                                   nComps,
                                   maxNeighbours,
                                   outputDir
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
                                     metaData = metaDataRef,
                                     samplesTrainDefList = createdModelsMajority$samplesTrainDefList,
                                     classColumn = classColumn,
                                     nModels = nModels,
                                     testSamples = testSamples,
                                     maxNeighbours = maxNeighbours
  )

  if (nrow(result) == 1) {
    randomVector <- paste0("fake", 1:ncol(result)) %>% as.data.frame() %>% t() %>% as.data.frame()
    colnames(randomVector) <- colnames(result)
    result1 <- rbind(result, randomVector)
    probability <-  apply(result1, 1, table)
  } else {
  # Find out how often a certain tumor type prediction is made for a specific sample
  probability <- apply(result, 1, table)
}
  # Locate the position of the highest probability
  positions <- lapply(probability, which.max)
  positions <- unlist(positions)

  bestFit <- data.frame(predict = rep(NA, times = length(result$fold1)))
  probabilityScores <- vector()

  # Extract the different calls being made for each sample
  mostAppearingNames <- lapply(probability, names)

  # Store the one with the highest probability score into the bestFit dataframe
  for (j in seq(1:length(mostAppearingNames))) {
    numberPositions <- as.numeric(positions[j])
    probabilityScores[j] <- probability[[j]][numberPositions]
    bestFit[j,] <- mostAppearingNames[[j]][numberPositions]
  }

  # Store the bestFit, the Originalcall and the accompanying probability score within the final dataframe.
  classifications <- cbind(bestFit,
                           probability = probabilityScores)

if (nrow(result) == 1) {
  classifications <- classifications[1, , drop = F]

}
  # Make sure that the classifications still have their accompanying biomaterial_id
  rownames(classifications) <- rownames(result)

  classificationList <- list(classifications = classifications,
                             probability = probability)

  directory <- paste0(outputDir, format(as.Date(Sys.Date(), "%Y-%m-%d"), "%m_%d_%Y"))
  if (!dir.exists(directory)) {
    dir.create(directory) }

  filename <- paste0(directory, "/majorityClassifierResult.rds")
  write_rds(classificationList, file = filename)

  return(classificationList)
}
