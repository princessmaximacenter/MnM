#' 10x cross-validation Minority classifier
#'
#' Setup to automatically run the 10x stratified cross-validation for the Minority classifier.
#' @param countDataRef Matrix containing the RNA-transcript per million data. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param metaDataRef Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param domainColumn Column in the metadata file that contains the tumor domain labels.
#' @param sampleColumn Column in the metadata file that contains the samples.
#' @param meanExpression Selection criterion for the RNA-transcripts,
#' specifying what the minimum mean expression of an RNA-transcript should be for it to be included in the F-statistic analysis.
#' Default is 5.
#' @param nANOVAgenes How many RNA-transcripts should we select using the F-statistic of ANOVA? Default is 1000.
#' @param maxSamplesPerType How many samples should we maximally use per tumor subtype? Default is 3.
#' @param nModels How many models should be created for the Minority classifier? Default is 100.
#' @param ntree How many trees should we use during the weighted Random Forest procedure? Default is 500.
#' @param nFeatures How many features should we keep after determining the most important RNA-transcripts using the Random Forest Importance Score?
#' Default is 300.
#' @param whichSeed For reproducibility, the seed can be specified with this parameter. Default is 1.
#' @param outputDir Directory in which you would like to store the R-object containing the results. Default is today's date.
#' @param proteinCodingGenes What are the names of the RNA-transcripts that stand for protein-coding genes within our dataset?
#' Please supply it as a vector. This is needed for ribo-depletion correction model.
#'
#' @return R-object containing the predictions ($classifications), classifications errors ($wrongClassifications),
#'  the probabilities for each classification ($probabilityList),
#'  the features that were eventually used for the weighted RF within the different folds ($reducedFeaturesList),
#'  the metadata file associated to the reference cohort ($metaDataRef)
#'  and the metadata for the performed Minority classifier run ($metaDataRun).
#' @export
#'@import foreach doParallel

tenFoldCrossValidationMinority <-  function(countDataRef,
                                            metaDataRef,
                                            classColumn,
                                            higherClassColumn,
                                            domainColumn,
                                            sampleColumn,
                                            meanExpression = 5,
                                            nANOVAgenes = 1000,
                                            maxSamplesPerType = 3,
                                            nModels = 100,
                                            ntree = 500,
                                            nFeatures = 300,
                                            whichSeed = 1,
                                            outputDir = paste0("./", format(as.Date(Sys.Date(), "%Y-%m-%d"), "%Y_%m_%d")),
                                            proteinCodingGenes

) {

  library(magrittr)
  `%notin%` <- Negate(`%in%`)

  if (sampleColumn %notin% colnames(metaDataRef)) {
    stop("The column you specified for the sample IDs is not present within metaDataRef. Please check the sampleColumn.")
  } else if (classColumn %notin% colnames(metaDataRef)) {
    stop("The column you specified for the tumor subtype labels is not present within metaDataRef. Please check the classColumn")
  } else if (higherClassColumn %notin% colnames(metaDataRef)){
    stop("The column you specified for the tumor type labels is not present within metaDataRef. Please check the higherClassColumn")
  } else if (domainColumn %notin% colnames(metaDataRef)) {
    stop("The column you specified for the tumor domain labels is not present within metaDataRef. Please check the domainColumn")
  }
  rownames(metaDataRef) <- metaDataRef[, sampleColumn]
  # Make sure the metadata and count data are in the right format and same order
  if (nrow(metaDataRef) != ncol(countDataRef)) {
    stop("The number of samples do not match between the metadata and the count data. Please make sure you include all same samples in both objects.")
  } else if (all(rownames(metaDataRef) %notin% colnames(countDataRef))) {
    stop("Your input data is not as required. Please make sure your patient IDs are within the row names of the metadata, and in the column names of the count data")
  }

  if (is.numeric(countDataRef) != T) {
    stop("Your input data is not as required. Please make sure your countDataRef object only contains numerical count data and is a matrix.
         Non-available measurements are not allowed.")

  }

  # Include a statement to store the classColumn, higherClassColumn and domainColumn
  print(paste0("The column used for tumor subtypes labels within the metadata, used for model training purposes, is: ", classColumn, ', containing values such as: '))
  print(base::unique(metaDataRef[,classColumn])[1:3])

  print(paste0("The column used for tumor type labels within the metadata, is: ", higherClassColumn,', containing values such as: '))
  print(base::unique(metaDataRef[,higherClassColumn])[1:3])

  print(paste0("The column used for tumor domain labels within the metadata, is: ", domainColumn, ', containing values such as: '))
  print(base::unique(metaDataRef[,domainColumn])[1:3])
  print("If any of these are incorrect, specify a different 'classColumn' (subtype), 'higherClassColumn' (tumor type) or 'domainColumn' (domain) to function as labels.")



  tumorEntitiesWithTooFewSamples <- base::table(metaDataRef[,classColumn])[base::table(metaDataRef[,classColumn]) < 3] %>% base::names()
  if (base::length(tumorEntitiesWithTooFewSamples) >0) {

    metaDataRef %<>% dplyr::filter(!!sym(classColumn) %notin% tumorEntitiesWithTooFewSamples)
    print("You have labels within your dataset that have less than 3 available samples. Please note samples with these labels have been removed.")

  }
  countDataRef <- countDataRef[, rownames(metaDataRef)]

  # Make sure you have CPM counts
  countDataRef <- apply(countDataRef,2,function(x) (x/sum(x))*1E6)

  if (!dir.exists(outputDir)) {
    dir.create(outputDir)}
  directory <- outputDir#paste0(outputDir, format(as.Date(Sys.Date(), "%Y-%m-%d"), "%m_%d_%Y"), "/")
  modelDirectory <- paste0(directory, "/seed", whichSeed)
  if (!dir.exists(directory)) {
    dir.create(directory)
    dir.create(modelDirectory)
  } else if (!dir.exists(modelDirectory)){
    dir.create(modelDirectory)
    }

  # Correct for ribosomal protein contamination
    riboCountFile <- paste0(modelDirectory, "/modelListRiboCounts.rds")
    if (!file.exists(riboCountFile)) {

      set.seed(whichSeed)
      riboModelList <- riboCorrectCounts(data = countDataRef,
                                         proteinCodingGenes = proteinCodingGenes,
                                         outputDir = modelDirectory,
                                         saveRiboModels = F
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
  folds <- caret::createFolds(metaDataRef[ , classColumn], k = 10, returnTrain = TRUE, list = TRUE)

  featuresAndModels <- foreach::foreach(i=c(1:base::length(folds))) %dopar%{

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
                                      nFeatures = nFeatures,
                                      nANOVAgenes = nANOVAgenes
                                      )

    dataCV <- dataCV[,c(reducedFeatures, "class")]

    #dataSynthList <- runUpsimpler("something", reducedFeatures = reducedFeatures)
    #dataSynth <- dataSynthList$dataSynth
    #metaDataSynth <- dataSynthList$metaDataSynth

    # Reduce features of testing data
    testDataCVReduced <- testDataCV[,reducedFeatures]
    print("Initiating RF")

    # Start the modelling of the data within the different compositions of training data
    set.seed(whichSeed)
    modelList <- obtainModelsMinorityClassifier(dataTrain = dataCV,
                                                #dataSynth = data
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

  for (i in seq(1:base::length(featuresAndModels))) {
    result <- featuresAndModels[[i]][["result"]]

    # Store (if required) models per fold within a list
    #modelList[[i]] <- featuresAndModels[[i]][["modelList"]]

    # Store reduced features per fold within a list
    reducedFeaturesList[[i]] <- featuresAndModels[[i]][["reducedFeatures"]]

    classificationResults <- convertResultToClassification(result = result,
                                                           metaDataRef = metaDataRef,
                                                           addOriginalCall = T,
                                                           classColumn = classColumn)

    probabilityList[[i]] <- classificationResults$probabilityList
    ultimatePredictions <- classificationResults$classifications

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
  accuracy <- sum(classifications$predict == classifications$originalCall) / base::length(classifications$originalCall)
  print(accuracy)

  # Store the settings of the classifier run within the resulting object
  metaDataRun <- data.frame(nModels = nModels,
                            classColumn = classColumn,
                            higherClassColumn = higherClassColumn,
                            domainColumn = domainColumn,
                            maxSamplesPerType = maxSamplesPerType,
                            nANOVAgenes = nANOVAgenes,
                            nFeatures = nFeatures,
                            whichSeed = whichSeed)

  # Save everything in the variable completePackage
  crossValidationMinorityResults <- list(classifications = classifications,
                                         wrongClassifications = wrongClassifications,
                                         probabilityList = probabilityList,
                                         reducedFeaturesList = reducedFeaturesList,
                                         metaDataRef = metaDataRef,
                                         metaDataRun = metaDataRun)

  print("We have finished the classification process. Please find your results in the generated object.")
  filename <- paste0(modelDirectory, "/crossValidationMinorityResults.rds")
  saveRDS(crossValidationMinorityResults, file = filename)
  print(paste0("Please find the generated R-object with the classification results within ", filename))
  return(crossValidationMinorityResults)
}
