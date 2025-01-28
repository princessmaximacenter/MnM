#' 10x cross-validation Minority classifier
#'
#' Setup to automatically run the 10x stratified cross-validation for the Minority classifier.
#' @param countDataRef Matrix containing the RNA-transcript per million data. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param metaDataRef Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Name of column in the metadata file that contains the tumor type labels.
#' @param domainColumn Name of column in the metadata file that contains the tumor domain labels.
#' @param sampleColumn Name of column in the metadata file that contains the sample identifiers.
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
#'@import foreach

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
                                            proteinCodingGenes,
                                            correctRibo = T

) {

  `%notin%` <- base::Negate(`%in%`)

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
  dataLogRef <- base::log(countDataRef +1) %>% base::t() %>% base::as.data.frame()

  # Specify max samples per tumor type
  if (base::is.na(maxSamplesPerType)) {
    maxSamplesPerType <- base::ceiling(stats::median(base::table(metaDataRef[, classColumn])))
  }

  # Remove genes with mean expression for ANOVA
  meanVals <- base::apply(countDataRef, 1, mean)
  countDataRef <- countDataRef[meanVals >= meanExpression,]

  cat("\nWe will now start with the cross-validation\n")
  # Set seed for reproducibility
  base::set.seed(whichSeed)

  # Create splits for cross-validation setup with equal distribution of tumor types
  folds <- caret::createFolds(metaDataRef[ , classColumn], k = 10, returnTrain = TRUE, list = TRUE)

  featuresAndModels <- foreach::foreach(i=c(1:base::length(folds))) %dopar%{

    # Select counts and metadata for training samples
    metaDataCV <- metaDataRef[folds[[i]],]
    countDataCV <- countDataRef[, folds[[i]]]

    # Run an ANOVA to select the top n genes from the training data for use in the further classification process
    base::set.seed(whichSeed)
    interestingANOVAgenes <- selectAnovaGenes(metaDataRef  = metaDataCV,
                                           countDataRef  = countDataCV,
                                           nANOVAgenes = nANOVAgenes,
                                           classColumn = classColumn)

    base::cat(base::paste("Found interesting ANOVA genes for Fold ", i, "\n"))
    # Select the ANOVA genes within the log-transformed data
    dataLogCV <- dataLogRef[ ,interestingANOVAgenes]

    # Select biomaterial IDs as training data per model
    base::set.seed(whichSeed)
    samplesTrainDefList <- obtainTrainData(metaDataRef = metaDataCV,
                                           classColumn = classColumn,
                                           nModels = nModels,
                                           maxSamplesPerType = maxSamplesPerType)
    # Select training data
    dataCV <- dataLogCV[folds[[i]], ]

    # Add class label to the dataset
    dataCV$class <- base::as.character(metaDataRef[rownames(dataCV),classColumn])

    # Also create test data and specify which biomaterial IDs are in there
    testDataCV <- dataLogCV[ -folds[[i]], ]
    #testSamples <- rownames(testDataCV)

    # Reduce features using RF feature importance for accuracy
    base::cat(base::paste("Working on reducing Features for Fold", i, "\n"))
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
    base::cat("\nInitiating RF\n")

    # Start the modelling of the data within the different compositions of training data
    base::set.seed(whichSeed)
    modelList <- obtainModelsMinorityClassifier(dataTrain = dataCV,
                                                #dataSynth = data
                                      samplesTrainDefList = samplesTrainDefList,
                                      nModels = nModels,
                                      ntree = ntree)

    base::cat(base::paste("Generated all models for Fold ", i, "\n"))
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
  base::cat("\nFinished with generating models and results\n")
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
  accuracy <- base::sum(classifications$predict == classifications$originalCall) / base::length(classifications$originalCall)
  base::cat(paste("\naccuracy: ", accuracy, "\n"))

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

  base::cat("\nWe have finished the classification process. Please find your results in the generated object.\n")
  filename <- base::paste0(modelDirectory, "/crossValidationMinorityResults.rds")
  base::saveRDS(crossValidationMinorityResults, file = filename)
  base::cat(base::paste0("Please find the generated R-object with the classification results within ", filename, "\n"))
  return(crossValidationMinorityResults)
}
