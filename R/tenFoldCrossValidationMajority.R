#' 10x cross-validation majority classifier
#'
#' Setup to automatically run the 10x stratified cross-validation majority classifier.
#'
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data. Patients are in the columns,
#' different genes in the rows.
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param nModels How many models should be created for the majority voting system?
#' @param maxSamplesPerType How many samples should we maximally use per tumor (sub)type?
#' @param nFeatures How many of the most variable genes within the dataset should we select for principal component analysis (PCA)?
#' @param nComps How many principal components will be selected after PCA?
#' @param maxNeighbours What is the maximum number of neigbours to be used for the weighted _k_-nearest neighbor algorithm?
#' @param whichSeed For reproducibility, the seed can be specified with this parameter.
#' @param outputDir Directory in which you would like to store the R-object containing the results.
#' @param proteinDir In which directory can we find the file specifying the names of protein-coding genes within our dataset?
#' @import tidyverse dplyr magrittr foreach doParallel kknn
#'
#' @return R-object containing the predictions ($classifications), classifications errors ($wrongClassifications),
#'  the probabilities for each classification ($probabilityList), the metadata file associated to the reference cohort ($metaData),
#'  and metadata for the performed run ($metaDataRun).
#' @export
#'
tenFoldCrossValidationMajority <-  function(countDataRef,
                                            metaDataRef,
                                            classColumn,
                                            nModels = 100,
                                            maxSamplesPerType = 50,
                                            nFeatures = 2500,
                                            nComps = 100,
                                            maxNeighbours = 25,
                                            whichSeed = 1,
                                            outputDir = "./",
                                            proteinDir = "~/surfdrive/Shared/Kemmeren group/Research_Projects/RNA_classification_FW/data/input/"

) {

  # Make sure you have CPM counts
  countDataRef <- apply(countDataRef,2,function(x) (x/sum(x))*1E6)


  directory <- paste0(outputDir, format(as.Date(Sys.Date(), "%Y-%m-%d"), "%m_%d_%Y"), "/")
  if (!dir.exists(directory)) {
    dir.create(directory) }

  # Correct for ribosomal protein contamination
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
    samplesTrainDefList <- obtainTrainData(metaDataRef = metaDataCV,
                                           classColumn = classColumn,
                                           nModels = nModels,
                                           maxSamplesPerType = maxSamplesPerType)


    rotationsAndScalingsList <- getPrincipalComponents(dataTrain = dataCV,
                                                       samplesTrainDefList,
                                                       classColumn = classColumn,
                                                       nModels = nModels,
                                                       nFeatures = nFeatures,
                                                       nComps = nComps
    )

    result <- obtainPredictionMajorityClassifier(rotationsAndScalingsList = rotationsAndScalingsList,
                                       dataTrain = dataCV,
                                       dataTest = testDataCV,
                                       metaDataRef = metaDataCV,
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
    originalCall <- metaDataRef[rownames(result),classColumn]

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
                                         metaData = metaDataRef,
                                         metaDataRun = metaDataRun
  )

  print("We have finished the classification process. Please find your results in the generated object.")
  filename <- paste0(directory, "/crossValidationMajorityResults.rds")
  write_rds(crossValidationMajorityResults, file = filename)
  return(crossValidationMajorityResults)
}
