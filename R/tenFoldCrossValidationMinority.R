#' 10x cross-validation minority classifier
#'
#' Setup to automatically run the 10x stratified cross-validation minority classifier.
#' @param countDataRef Matrix containing the RNA-transcript per million data. Patients are in the columns,
#' different genes in the rows.
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param meanExpression Selection criterion for the genes,
#' specifying what the minimum mean expression of a gene should be for it to be included in the F-statistic analysis.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param nModels How many models should be created for the majority voting system?
#' @param nANOVAgenes How many genes should we select using the F-statistic of ANOVA?
#' @param maxSamplesPerType How many samples should we maximally use per tumor (sub)type?
#' @param ntree How many trees should we use during the weighted Random Forest (RF) procedure?
#' @param howManyFeatures How many features should we keep after determining the most important genes using the Random Forest Importance Score?
#' @param whichSeed For reproducibility, the seed can be specified with this parameter.
#' @param outputDir Directory in which you would like to store the R-object containing the results.
#' @param proteinDir In which directory can we find the file specifying the names of protein-coding genes within our dataset?
#'
#' @return R-object containing the predictions ($classifications), classifications errors ($wrongClassifications),
#'  the probabilities for each classification ($probabilityList), the metadata file associated to the reference cohort ($metaData),
#' the features that were eventually used for the weighted RF within the different folds ($reducedFeaturesList)
#'  and the metadata for the performed run ($metaDataRun).
#' @export
#'@import tidyverse dplyr magrittr foreach doParallel randomForest caret glmnet

tenFoldCrossValidationMinority <-  function(countDataRef,
                                            metaDataRef,
                                            meanExpression = 5,
                                            classColumn,
                                            nModels = 100,
                                            nANOVAgenes = 1000,
                                            maxSamplesPerType = 3,
                                            ntree = 500,
                                            howManyFeatures = 300,
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
  dataLogRef <- log(countDataRef +1) %>% t() %>% as.data.frame()

  # Specify max samples per tumor type
  if (is.na(maxSamplesPerType)) {
    maxSamplesPerType <- ceiling(median(table(metaDataRef[, classColumn])))
  }

  # Remove genes with mean expression for ANOVA
  meanVals <- apply(countDataRef, 1, mean)
  countDataRef <- countDataRef[meanVals >= meanExpression,]

  print("We will now start with the cross-validation")
  # Set seed for reproducibility
  set.seed(whichSeed)

  # Create splits for cross-validation setup with equal distribution of tumor types
  folds <- createFolds(metaDataRef[ , classColumn], k = 10, returnTrain = TRUE, list = TRUE)

  featuresAndModels <- foreach(i=c(1:length(folds))) %dopar%{

    # Select counts and metadata for training samples
    metaDataCV <- metaDataRef[folds[[i]],]
    countDataCV <- countDataRef[, folds[[i]]]

    # Run an ANOVA to select the top n genes from the training data for use in the further classification process
    set.seed(whichSeed)
    interestingANOVAgenes <- selectAnovaGenes(metaDataRef  = metaDataCV,
                                           countDataRef  = countDataCV,
                                           nANOVAgenes = nANOVAgenes, # How many ANOVA genes
                                           classColumn = classColumn)

    print(paste("Found interesting ANOVA genes for Fold ", i))
    # Select the ANOVA genes within the log-transformed data
    dataLogCV <- dataLogRef[ ,interestingANOVAgenes]

    # Select biomaterial IDs as training data per model
    set.seed(whichSeed)
    samplesTrainDefList <- obtainTrainData(metaDataRef = metaDataCV,
                                           classColumn = classColumn,
                                           nModels = nModels,
                                           maxSamplesPerType = maxSamplesPerType)
    # Select training data
    dataCV <- dataLogCV[folds[[i]], ]

    # Add class label to the dataset
    dataCV$class <- as.character(metaDataRef[rownames(dataCV),classColumn])

    # Also create test data and specify which biomaterial IDs are in there
    testDataCV <- dataLogCV[ -folds[[i]], ]
    #testSamples <- rownames(testDataCV)

    # Reduce features using RF feature importance for accuracy
    print(paste("Working on reducing Features for Fold", i))
    reducedFeatures <- reduceFeatures(dataTrain = dataCV,
                                      samplesTrainDefList = samplesTrainDefList,
                                      ntree = 500,
                                      nModels = nModels,
                                      howManyFeatures = howManyFeatures)

    dataCV <- dataCV[,c(reducedFeatures, "class")]

    # Reduce features of testing data
    testDataCVReduced <- testDataCV[,reducedFeatures]
    print("Initiating RF")

    # Start the modelling of the data within the different compositions of training data
    set.seed(whichSeed)
    modelList <- obtainModelsMinorityClassifier(dataTrain = dataCV,
                                      samplesTrainDefList = samplesTrainDefList,
                                      nModels = nModels,
                                      ntree = ntree)

    print(paste("Generated all models for Fold ", i))
    # Find the predictions for the test data
    result <- predictTest(modelList, testDataCVReduced)

    # Combine results, generated models and selected features into a list
    featuresAndModels <- list(result = result,
                              #modelList = modelList,
                              reducedFeatures = reducedFeatures)

    # In the end you will end up with the resulting classifications per tumor sample,
    # the generated models and if feature reduction took place the features chosen per fold.
    return(featuresAndModels)
  }
  print("Finished with generating models and results")
  probabilityList <- list()

  #modelList <- list()
  reducedFeaturesList <- list()

  # We need to combine the results for the different folds.
  # This is done looping over all folds.

  for (i in seq(1:length(featuresAndModels))) {
    result <- featuresAndModels[[i]][["result"]]

    # Store (if required) models per fold within a list
    #modelList[[i]] <- featuresAndModels[[i]][["modelList"]]

    # Store reduced features per fold within a list
    reducedFeaturesList[[i]] <- featuresAndModels[[i]][["reducedFeatures"]]

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
                            nANOVAgenes = nANOVAgenes,
                            howManyFeatures = howManyFeatures,
                            whichSeed = whichSeed)

  # Save everything in the variable completePackage
  crossValidationMinorityResults <- list(classifications = classifications,
                                         wrongClassifications = wrongClassifications,
                                         probabilityList = probabilityList,
                                         reducedFeaturesList = reducedFeaturesList,
                                         metaData = metaDataRef,
                                         metaDataRun = metaDataRun)

  print("We have finished the classification process. Please find your results in the generated object.")
  filename <- paste0(directory, "/crossValidationMinorityResults.rds")
  write_rds(crossValidationMinorityResults, file = filename)
  return(crossValidationMinorityResults)
}
