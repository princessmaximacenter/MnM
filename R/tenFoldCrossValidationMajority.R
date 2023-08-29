tenFoldCrossValidationMajority <-  function(countDataRef,
                                            metaDataRef,
                                            classColumn,
                                            nModels = 100,
                                            maxSamplesPerType = 50,
                                            nComps = 100,
                                            nFeatures = 2500,
                                            maxNeighbours = 25,
                                            whichSeed = 1,
                                            outputDir = "./",
                                            proteinDir = "~/surfdrive/Shared/Kemmeren group/Research_Projects/RNA_classification_FW/data/input/",
                                            correctRibosomes = F

) {

  # Make sure you have CPM counts
  countDataRef <- apply(countDataRef,2,function(x) (x/sum(x))*1E6)


  directory <- paste0(outputDir, format(as.Date(Sys.Date(), "%Y-%m-%d"), "%m_%d_%Y"), "/")
  if (!dir.exists(directory)) {
    dir.create(directory) }

  # Correct for ribosomal protein contamination
  if (correctRibosomes == T) {
    riboCountFile <- paste0(directory, "modelListRiboCounts.rds")
    if (!file.exists(riboCountFile)) {

      proteinCodingGenes <- read.table(paste0(proteinDir,"20230320_proteinCodingGenes_gencode31.csv"), sep = "\t") %>%
        select(x) %>% deframe
      set.seed(whichSeed)
      riboModelList <- riboCorrectCounts(data = countDataRef,
                                         proteinCodingGenes = proteinCodingGenes,
                                         outputDir = directory
      )

    } else {
      riboModelList <- readRDS(riboCountFile)

    }
    countDataRef <- riboModelList$counts
  }

  # Log-transform data
  dataLogRef <- log(countDataRef +1)

  # Remove the genes with zero variance across the dataset

  dataLogZeroVar <- t(dataLogRef) %>% as.data.frame(.)
  zeroVar <- nearZeroVar(dataLogZeroVar)

  dataLogNonZero <- dataLogRef[-zeroVar, ]

  print("We will now start with the cross-validation")
  # Set seed for reproducibility
  set.seed(whichSeed)

  # Create splits for cross-validation setup with equal distribution of tumor types
  folds <- createFolds(metaDataRef[ , classColumn], k = 10, returnTrain = TRUE, list = TRUE)

  featuresAndModels <- foreach(i=c(1:length(folds))) %dopar%{

    print(paste0("Working on fold", i))
    # Select metadata and log-transformed data for training samples
    metaDataCV <- metaDataRef[folds[[i]],]
    dataCV <- dataLogNonZero[ , folds[[i]]]


    # Also create test data and specify which biomaterial IDs are in there
    testDataCV <- dataLogNonZero[ , -folds[[i]]]
    testSamples <- colnames(testDataCV)

    # Select biomaterial IDs as training data per model
    set.seed(whichSeed)
    samplesTrainDefList <- ObtainTrainData(metaData = metaDataCV,
                                           classColumn = classColumn,
                                           nModels = nModels,
                                           maxSamplesPerType = maxSamplesPerType)


    rotationsAndScalingsList <- getPrincipalComponents(dataTrain = dataCV,
                                                       metaData = metaDataRef,
                                                       samplesTrainDefList,
                                                       classColumn = classColumn,
                                                       nModels = nModels,
                                                       nFeatures = nFeatures,
                                                       nComps = nComps
    )

    result <- obtainPredictionMajority(rotationsAndScalingsList = rotationsAndScalingsList,
                                       dataTrain = dataCV,
                                       dataTest = testDataCV,
                                       metaData = metaDataRef,
                                       samplesTrainDefList = samplesTrainDefList,
                                       testSamples = testSamples,
                                       classColumn = classColumn,
                                       nModels = nModels,
                                       nComps = nComps,
                                       maxNeighbours = maxNeighbours
    )

    featuresAndModels <- list(result = result)

    return(featuresAndModels)
  }

  probabilityList <- list()

  # We need to combine the results for the different folds.
  # This is done looping over all folds.

  for (i in seq(1:length(featuresAndModels))) {
    result <- featuresAndModels[[i]][["result"]]

    # Find out how often a certain tumor type prediction is made for a specific sample
    probability <- apply(result, 1, table)

    # Store the resulting probabilities in a list
    probabilityList[[i]] <- probability

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

    # Look at the original calls for each test sample
    originalCall <- metaData[rownames(result),classColumn]

    # Store the bestFit, the Originalcall and the accompanying probability score within the final dataframe.
    ultimatePredictions <- cbind(bestFit,
                                 originalCall = originalCall,
                                 probability = probabilityScores)


    # Make sure that the ultimatePredictions still have their accompanying biomaterial_id
    rownames(ultimatePredictions) <- rownames(result)

    # Add each fold to a dataframe that combines all results
    if (i == 1) {
      wrongClassifications <- ultimatePredictions[ultimatePredictions$predict != ultimatePredictions$originalCall,]
      classifications <- ultimatePredictions
    } else {
      wrongClassifications <- rbind(wrongClassifications, ultimatePredictions[ultimatePredictions$predict != ultimatePredictions$originalCall,])
      classifications <- rbind(classifications, ultimatePredictions)
    }
  }

  # Check the accuracy of the current run
  accuracy <- sum(classifications$predict == classifications$originalCall) / length(classifications$originalCall)
  print(accuracy)

  # Store the settings of the classifier run within the resulting object
  metaDataRun <- data.frame(nModels = nModels,
                            classColumn = classColumn,
                            maxSamplesPerType = maxSamplesPerType,
                            nFeatures = nFeatures,
                            nComps = nComps,
                            whichSeed = whichSeed)

  # Save everything in the variable completePackage
  crossValidationMajorityResults <- list(classifications = classifications,
                                         wrongClassifications = wrongClassifications,
                                         probabilityList = probabilityList,
                                         metaData = metaData,
                                         metaDataRun = metaDataRun
  )

  print("We have finished the classification process. Please find your results in the generated object.")
  filename <- paste0(directory, "/crossValidationMajorityResults.rds")
  write_rds(crossValidationMajorityResults, file = filename)
  return(crossValidationMajorityResults)
}
