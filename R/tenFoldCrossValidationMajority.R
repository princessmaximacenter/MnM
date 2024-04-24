#' 10x cross-validation majority classifier
#'
#' Setup to automatically run the 10x stratified cross-validation majority classifier.
#'
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param metaDataRef Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param domainColumn Column in the metadata file that contains the tumor domain labels.
#' @param sampleColumn Column in the metadata file that contains the samples.
#' @param nModels How many models should be created for the Majority classifier?
#' @param maxSamplesPerType How many samples should we maximally use per tumor (sub)type?
#' @param nFeatures How many of the most variable genes within the dataset should we select for principal component analysis (PCA)?
#' @param nComps How many principal components will be selected after PCA?
#' @param maxNeighbors What is the maximum number of neigbours to be used for the weighted _k_-nearest neighbor algorithm?
#' @param whichSeed For reproducibility, the seed can be specified with this parameter.
#' @param outputDir Directory in which you would like to store the R-object containing the results. Default is today's date.
#' @param proteinCodingGenes What are the names of the RNA-transcripts that stand for protein-coding genes within our dataset?
#' Please supply it as a vector. This is needed for ribo-depletion correction model.
#' @import foreach doParallel
#'
#' @return R-object containing the predictions ($classifications), classifications errors ($wrongClassifications),
#'  the probabilities for each classification ($probabilityList), the  rotations and scalings
#' for each reference cohort subset($rotationsAndScalingList),
#' the metadata file associated to the reference cohort ($metaDataRef),
#'  and metadata for the performed run ($metaDataRun).
#' @export
#'
tenFoldCrossValidationMajority <-  function(countDataRef,
                                            metaDataRef,
                                            classColumn,
                                            higherClassColumn,
                                            domainColumn,
                                            sampleColumn,
                                            nModels = 100,
                                            maxSamplesPerType = 50,
                                            nFeatures = 2500,
                                            nComps = 100,
                                            maxNeighbors = 25,
                                            whichSeed = 1,
                                            outputDir = paste0("./", format(as.Date(Sys.Date(), "%Y-%m-%d"), "%Y_%m_%d")),
                                            proteinCodingGenes

) {


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
    stop("Your input data is not as required. Please make sure your sample IDs are within the sampleColumn or within the row names of the metadata, and in the column names of the count data")
  }

  if (is.numeric(countDataRef) != T) {
    stop("Your input data is not as required. Please make sure your countDataRef object only contains numerical count data and is a matrix.")

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
  if (length(tumorEntitiesWithTooFewSamples) >0) {

    metaDataRef %<>% filter(!!sym(classColumn) %notin% tumorEntitiesWithTooFewSamples)
    print("You have labels within your dataset that have less than 3 available samples.  Please note samples with these labels have been removed.")
    #stop("You have tumor subtypes within your dataset that have less than 3 available samples. Please remove all tumor types with too few samples. ")

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
  dataLogRef <- log(countDataRef +1)

  # Remove the genes with zero variance across the dataset

  dataLogZeroVar <- t(dataLogRef) %>% as.data.frame(.)
  zeroVar <- caret::nearZeroVar(dataLogZeroVar)

  dataLogNonZero <- dataLogRef[-zeroVar, ]

  print("We will now start with the cross-validation")
  # Set seed for reproducibility
  set.seed(whichSeed)

  # Create splits for cross-validation setup with equal distribution of tumor types
  folds <- caret::createFolds(metaDataRef[ , classColumn], k = 10, returnTrain = TRUE, list = TRUE)

  featuresAndModels <- foreach::foreach(i=c(1:base::length(folds))) %dopar%{

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
                                       maxNeighbors = maxNeighbors
    )

    featuresAndModels <- list(result = result)

    return(featuresAndModels)
  }

  probabilityList <- list()

  # We need to combine the results for the different folds.
  # This is done looping over all folds.

  for (i in seq(1:base::length(featuresAndModels))) {
    result <- featuresAndModels[[i]][["result"]]


    #FUNCTION TO GET THE ultimatePredictions"
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
                            maxNeighbors = maxNeighbors,
                            nFeatures = nFeatures,
                            nComps = nComps,
                            whichSeed = whichSeed)

  # Save everything in the variable completePackage
  crossValidationMajorityResults <- list(classifications = classifications,
                                         wrongClassifications = wrongClassifications,
                                         probabilityList = probabilityList,
                                         rotationsAndScalingsList = rotationsAndScalingsList,
                                         metaDataRef = metaDataRef,
                                         metaDataRun = metaDataRun
  )

  print("We have finished the classification process. Please find your results in the generated object.")
  filename <- paste0(modelDirectory, "/crossValidationMajorityResults.rds")
  saveRDS(crossValidationMajorityResults, file = filename)
  print(paste0("Please find the generated R-object with the classification results within ", filename))
  return(crossValidationMajorityResults)
}
