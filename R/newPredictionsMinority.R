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


  if (nrow(result) == 1 || typeof(apply(result, 1, table)) == 'integer') {
    randomVector <- paste0("fake", 1:ncol(result)) %>% as.data.frame() %>% t() %>% as.data.frame()
    colnames(randomVector) <- colnames(result)
    result1 <- rbind(result, randomVector)
    probability <-  apply(result1, 1, table)
    probability <- probability[-length(probability)]
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
                             probabilityList = probability,
                             metaDataRef = createdModelsMinority$metaDataRef)
  if (saveModel == T) {
    #directory <- outputDir
    filename <- paste0(outputDir, "/minorityClassifierResult.rds")
    if (!dir.exists(outputDir)) {
      dir.create(outputDir) }
  saveRDS(classificationList, file = filename)
  }
  return(classificationList)
}
