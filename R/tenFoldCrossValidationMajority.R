#' 10x cross-validation majority classifier
#'
#' Setup to automatically run the 10x stratified cross-validation majority classifier.
#'
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param metaDataRef Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Name of column in the metadata file that contains the tumor type labels.
#' @param domainColumn Name of column in the metadata file that contains the tumor domain labels.
#' @param sampleColumn Name of column in the metadata file that contains the sample identifiers.
#' @param nModels How many models should be created for the Majority classifier?
#' @param maxSamplesPerType How many samples should we maximally use per tumor (sub)type?
#' @param nFeatures How many of the most variable genes within the dataset should we select for principal component analysis (PCA)?
#' @param nComps How many principal components will be selected after PCA?
#' @param maxNeighbors What is the maximum number of neighbors to be used for the weighted _k_-nearest neighbor algorithm?
#' @param whichSeed For reproducibility, the seed can be specified with this parameter.
#' @param outputDir Directory in which you would like to store the R-object containing the results. Default is today's date.
#' @param proteinCodingGenes What are the names of the RNA-transcripts that stand for protein-coding genes within our dataset?
#' Please supply it as a vector. This is needed for ribo-depletion correction model.
#' @import foreach
#'
#' @return R-object containing the predictions ($classifications), classifications errors ($wrongClassifications),
#'  the probabilities for each classification ($probabilityList),
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
                                            proteinCodingGenes,
                                            correctRibo = T

) {


  `%notin%` <<- base::Negate(`%in%`)

  checkFormatInputData(sampleColumn = sampleColumn,
                       classColumn = classColumn,
                       higherClassColumn = higherClassColumn,
                       domainColumn = domainColumn,
                       metaDataRef = metaDataRef,
                       countDataRef = countDataRef,
                       outputDir = outputDir,
                       saveModel = T)

  base::rownames(metaDataRef) <- metaDataRef[, sampleColumn]

  tumorEntitiesWithTooFewSamples <- base::table(metaDataRef[,classColumn])[base::table(metaDataRef[,classColumn]) < 3] %>% base::names()

  if (base::length(tumorEntitiesWithTooFewSamples) >0) {

    metaDataRef %<>% dplyr::filter(!!dplyr::sym(classColumn) %notin% tumorEntitiesWithTooFewSamples)
    base::cat("\nYou have labels within your dataset that have less than 3 available samples. \nPlease note samples with these labels have been removed.\n")
    #stop("You have tumor subtypes within your dataset that have less than 3 available samples. Please remove all tumor types with too few samples. ")

  }
  countDataRef <- countDataRef[, base::rownames(metaDataRef)]

  # Make sure you have CPM counts
  countDataRef <- base::apply(countDataRef,2,function(x) (x/base::sum(x))*1E6)

  directory <- outputDir#paste0(outputDir, format(as.Date(Sys.Date(), "%Y-%m-%d"), "%m_%d_%Y"), "/")
  modelDirectory <- base::paste0(directory, "/seed", whichSeed)

  if (!base::dir.exists(modelDirectory)){
    base::dir.create(modelDirectory)
  }

  if(correctRibo == T) {
    # Correct for ribosomal protein contamination
    riboCountFile <- base::paste0(directory, "modelListRiboCounts.rds")
    if (!base::file.exists(riboCountFile)) {

      base::set.seed(whichSeed)
      riboModelList <- riboCorrectCounts(data = countDataRef,
                                         proteinCodingGenes = proteinCodingGenes,
                                         outputDir = directory
      )

    } else {
      riboModelList <- base::readRDS(riboCountFile)
    }
    base::cat("\nFinished the ribocorrection\n")
    countDataRef <- riboModelList$counts
  }

  # Log-transform data
  dataLogRef <- base::log(countDataRef +1)

  # Remove the genes with zero variance across the dataset

  dataLogZeroVar <- base::t(dataLogRef) %>% base::as.data.frame()
  zeroVar <- caret::nearZeroVar(dataLogZeroVar)

  dataLogNonZero <- dataLogRef[-zeroVar, ]

  base::print("We will now start with the cross-validation")
  # Set seed for reproducibility
  base::set.seed(whichSeed)

  # Create splits for cross-validation setup with equal distribution of tumor types
  folds <- caret::createFolds(metaDataRef[ , classColumn], k = 10, returnTrain = TRUE, list = TRUE)

  featuresAndModels <- foreach::foreach(i=c(1:base::length(folds))) %dopar%{

    base::print(base::paste0("Working on fold", i))
    # Select metadata and log-transformed data for training samples
    metaDataCV <- metaDataRef[folds[[i]],]
    dataCV <- dataLogNonZero[ , folds[[i]]]


    # Also create test data and specify which biomaterial IDs are in there
    testDataCV <- dataLogNonZero[ , -folds[[i]], drop = F]
    testSamples <- colnames(testDataCV)

    # Select biomaterial IDs as training data per model
    base::set.seed(whichSeed)
    samplesTrainDefList <- obtainTrainData(metaDataRef = metaDataCV,
                                           classColumn = classColumn,
                                           nModels = nModels,
                                           maxSamplesPerType = maxSamplesPerType)


    # This function needs to be split up into 2
    scaleFeaturesList <- getVarFeaturesMajority(dataTrain = dataLogNonZero,
                                                samplesTrainDefList = samplesTrainDefList,
                                                nFeatures = nFeatures,
                                                nModels = nModels )

    base::cat("\nSelected features and calculated scaling\n")


    rotationsAndScalingsList <- getPrincipalComponents(dataTrain = dataLogNonZero,
                                                       scaleFeaturesList = scaleFeaturesList,
                                                       samplesTrainDefList = samplesTrainDefList,
                                                       nModels = nModels,
                                                       nComps = nComps)

    vectorK <- trainKNNClassifier(rotationsAndScalingsList = rotationsAndScalingsList,
                                  metaDataRef = metaDataCV,
                                  classColumn = classColumn,
                                  nModels = nModels,
                                  maxNeighbors = maxNeighbors)

    result <- obtainPredictionMajorityClassifier(rotationsAndScalingsList = rotationsAndScalingsList,
                                       dataTest = testDataCV,
                                       metaDataRef = metaDataCV,
                                       testSamples = testSamples,
                                       classColumn = classColumn,
                                       nModels = nModels,
                                       vectorK = vectorK
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
  accuracy <- base::sum(classifications$predict == classifications$originalCall) / base::length(classifications$originalCall)
  base::print(base::paste("accuracy: ", accuracy))

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
                                         metaDataRef = metaDataRef,
                                         metaDataRun = metaDataRun
  )

  base::print("We have finished the classification process. Please find your results in the generated object.")
  filename <- base::paste0(modelDirectory, "/crossValidationMajorityResults.rds")
  base::saveRDS(crossValidationMajorityResults, file = filename)
  base::print(base::paste0("Please find the generated R-object with the classification results within ", filename))
  return(crossValidationMajorityResults)
}
