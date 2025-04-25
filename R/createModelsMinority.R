#' Generate Minority classifier models
#'
#' Setup to obtain the models for the Minority classifier.
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param metaDataRef Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Name of column in the metadata file that contains the tumor type labels.
#' @param domainColumn Name of column in the metadata file that contains the tumor domain labels.
#' @param sampleColumn Name of column in the metadata file that contains the sample identifiers.
#' @param nModels How many models should be created for the Minority classifier?
#' Default is 100.
#' @param meanExpression Selection criterion for the RNA-transcripts,
#' specifying what the minimum mean TPM-expression of a RNA-transcripts should be for it to be included in the F-statistic analysis.
#' Default is 5.
#' @param nANOVAgenes How many RNA-transcripts should we select using the F-statistic of ANOVA? Default is 1000.
#' @param maxSamplesPerType How many samples should we maximally use per tumor subtype? Default is 3.
#' @param ntree How many trees should we use during the weighted Random Forest procedure? Default is 500.
#' @param nFeatures How many features should we keep after determining the most important RNA-transcripts using the Random Forest Importance Score?
#' Default is 300.
#' @param whichSeed For reproducibility, the seed can be specified with this parameter. Default is 1.
#' @param outputDir Directory in which you would like to store the R-object containing the results. Default is today's date.
#' @param proteinCodingGenes What are the names of the RNA-transcripts that stand for protein-coding genes within our dataset?
#' Please supply it as a vector. This is needed for ribo-depletion correction model.
#' @param saveModel Do you want to save your generated scalings in an R-object? Default is TRUE.
#' @param useUpsimpler Do you want to use the Upsimpler modality to create new samples?
#' @param upsampleBelow Which is the number of samples you want to leave as is within M&M, below which upsampling will be performed?
#' @param upsamplerArgs What are the parameters you would like to use for upsampling?
#' @param upsamplerModule Module containing the upsampling package.
#' @param upsamplingType Do you want to use upsimpler or random upsampling?
#' @param featuresToRemove Do you want to remove specific features so they are not being considered for use within the classification procedure?
#' @return R-object containing the generated RF-models ($modelList), the model for the ribodepletion correction ($riboModelList),
#'  the features that were eventually used for the weighted RF within the different folds ($reducedFeaturesList),
#'   the metadata file associated to the reference cohort ($metaDataRef)
#'  and the metadata for the performed run ($metaDataRun).
#'
#' @export
#'
createModelsMinority <-  function(countDataRef,
                                  metaDataRef,
                                  classColumn,
                                  higherClassColumn,
                                  domainColumn,
                                  sampleColumn,
                                  nModels = 100,
                                  meanExpression = 5,
                                  nANOVAgenes = 1000,
                                  maxSamplesPerType = 3,
                                  ntree = 500,
                                  nFeatures = 300,
                                  whichSeed = 1,
                                  outputDir = paste0("./", format(as.Date(Sys.Date(), "%Y-%m-%d"), "%Y_%m_%d")),
                                  proteinCodingGenes,
                                  saveModel = T,
                                  correctRibo = T,
                                  useUpsimpler = F,
                                  upsampleBelow = 4,
                                  maxSamplesAfterUpsampling = 6,
                                  upsamplerArgs = NULL,
                                  upsamplerModule = NULL,
                                  upsamplingType = "upsimpler",
                                  featuresToRemove = NULL

) {
  countDataOG <- countDataRef
  `%notin%` <<- base::Negate(`%in%`)

  checkFormatInputData(sampleColumn = sampleColumn,
                       classColumn = classColumn,
                       higherClassColumn = higherClassColumn,
                       domainColumn = domainColumn,
                       metaDataRef = metaDataRef,
                       countDataRef = countDataRef,
                       outputDir = outputDir,
                       saveModel = saveModel)

  tumorEntitiesWithTooFewSamples <- base::table(metaDataRef[,classColumn])[base::table(metaDataRef[,classColumn]) < 3] %>% base::names()
  base::rownames(metaDataRef) <- metaDataRef[, sampleColumn]

  if (base::length(tumorEntitiesWithTooFewSamples) >0 & useUpsimpler == F) {

    metaDataRef %<>% dplyr::filter(!!dplyr::sym(classColumn) %notin% tumorEntitiesWithTooFewSamples)
    base::print("You have labels within your dataset that have less than 3 available samples.  Please note samples with these labels have been removed.")
    #stop("You have tumor subtypes within your dataset that have less than 3 available samples. Please remove all tumor types with too few samples. ")

  }

  countDataRef <- countDataRef[, base::rownames(metaDataRef)]
  # Make sure you have CPM counts
  countDataRef <- base::apply(countDataRef,2,function(x) (x/base::sum(x))*1E6)

  directory <- outputDir

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

  # Remove features if needed
  if (!is.null(featuresToRemove)) {
    # Remove features by rownames.
    countDataRef <- countDataRef %>% base::as.data.frame() %>%
      dplyr::filter(rownames(.) %notin% featuresToRemove) %>%
      base::as.matrix()
  }

  # Log-transform data
  dataLogRef <- base::log(countDataRef +1) %>% base::t() %>% base::as.data.frame()

  # Specify max samples per tumor type
  if (base::is.na(maxSamplesPerType)) {
    maxSamplesPerType <- base::ceiling(stats::median(base::table(metaDataRef[, classColumn])))
  }

  # Remove genes with mean expression for ANOVA
  meanVals <- base::apply(countDataRef, 1, base::mean)
  countDataRef <- countDataRef[meanVals >= meanExpression,]

  base::print("We will now start with the selection of ANOVA-genes")
  # Set seed for reproducibility
  base::set.seed(whichSeed)


  # Run an ANOVA to select the top n genes from the training data for use in the further classification process
  interestingANOVAgenes <- selectAnovaGenes(metaDataRef = metaDataRef,
                                         countDataRef = countDataRef,
                                         nANOVAgenes = nANOVAgenes, # How many ANOVA genes
                                         classColumn = classColumn)

  base::print("We have finished with the ANOVA gene selection")
  # Select the ANOVA genes within the log-transformed data
  dataLogRef <- dataLogRef[ ,interestingANOVAgenes]

  # Select biomaterial IDs as training data per model
  base::set.seed(whichSeed)
  samplesTrainDefList <- obtainTrainData(metaDataRef = metaDataRef,
                                         classColumn = classColumn,
                                         nModels = nModels,
                                         maxSamplesPerType = maxSamplesPerType)

  # Add class label to the dataset
  dataLogRef$class <- base::as.character(metaDataRef[rownames(dataLogRef),classColumn])

  base::print("Starting to reduce features")
  # Reduce features using RF feature importance for accuracy
  base::set.seed(whichSeed)
  reducedFeatures <- reduceFeatures(dataTrain = dataLogRef,
                                    samplesTrainDefList = samplesTrainDefList,
                                    ntree = 500,
                                    nModels = nModels,
                                    nFeatures = nFeatures,
                                    nANOVAgenes = nANOVAgenes)

  dataLogRef <- dataLogRef[,c(reducedFeatures, "class")]

  tumorEntitiesWithTooFewSamples <- base::table(metaDataRef[,classColumn])[base::table(metaDataRef[,classColumn]) < upsampleBelow] %>% base::names()
  # Apply upsampling to create new samples
  if (useUpsimpler == T & length(tumorEntitiesWithTooFewSamples) >= 1) {
    #tumorEntitiesWithTooFewSamples <- base::table(metaDataRef[,classColumn])[base::table(metaDataRef[,classColumn]) < upsampleBelow] %>% base::names()
    #cat(tumorEntitiesWithTooFewSamples)
    seedSamples <- metaDataRef %>% dplyr::filter(!!sym(classColumn) %in% tumorEntitiesWithTooFewSamples) %>%
      dplyr::pull(!!dplyr::sym(sampleColumn))

    if(is.null(upsamplerArgs)) {
      upsamplerArgs <- getDefaultSettingsUpsampling(upsamplingType = upsamplingType)
    }

    if (upsamplingType == "upsimpler") {
      upsamplerArgs$init$dataDF <- dataLogRef %>% dplyr::select(!class)
      upsamplerArgs$init$metadataDF <- metaDataRef
      upsamplerArgs$init$class_col <- classColumn
      upsimpler <- base::do.call(upsamplerModule$Upsimpler, upsamplerArgs$init)


    } else if (upsamplingType == "random") {
      upsamplerArgs$dataDF <- dataLogRef %>% dplyr::select(!class)
      upsamplerArgs$metadataDF <- metaDataRef
      upsamplerArgs$class_col <- classColumn
      upsimpler <- upsamplerModule$baseline

    }



      synths <- applyUpsimpler(upsimpler = upsimpler,
                               upsamplingType = upsamplingType,
                               upsimplerUsedArgs = upsamplerArgs,
                               targetSampleIDs = seedSamples,
                               metadataDF = metaDataRef,
                               classColumn = classColumn,
                               domainColumn = domainColumn,
                               higherClassColumn = higherClassColumn,
                               sampleColumn = sampleColumn)

      countsSynths <- base::expm1(base::t(synths$synthDataDF)) # For reference later on

      if (upsamplingType == "upsimpler") {
        dataLogRef <- base::rbind(upsamplerArgs$init$dataDF, synths$synthDataDF)
      } else if (upsamplingType == "random") {
        dataLogRef <- base::rbind(upsamplerArgs$dataDF, synths$synthDataDF)
      }

      # synths metadata is a single column; we need to add the other columns
      # in order to rbind it to the original metadata
      synths$synthMetadataDF[base::setdiff(base::names(metaDataRef), base::names(synths$synthMetadataDF))] <- NA
      metaDataWithSynths <- base::rbind(metaDataRef, synths$synthMetadataDF)

    # Alter training dataset
    samplesTrainDefList <- obtainTrainData(metaDataRef = metaDataWithSynths,
                    classColumn = classColumn,
                    nModels = nModels,
                    maxSamplesPerType = maxSamplesAfterUpsampling)
    # Update the samplesTrainDefList to reflect the new dataLogRef (i.e., including the synths),
    # by breaking the synths up into the different models
    # synthSamplesTrainDefList <- obtainTrainData(metaDataRef = synths$synthMetadataDF,
    #                                             classColumn = classColumn,
    #                                             nModels = nModels,
    #                                             maxSamplesPerType = maxSamplesPerType)
    # samplesTrainDefList <- base::Map(base::c, samplesTrainDefList, synthSamplesTrainDefList)

    # Add back the class label to the dataset
    dataLogRef$class <- base::as.character(metaDataWithSynths[rownames(dataLogRef), classColumn])
  }



  base::print("Initiating RF")

  # Start the modelling of the data within the different compositions of training data
  base::set.seed(whichSeed)
  modelList <- obtainModelsMinorityClassifier(dataTrain = dataLogRef,
                                    samplesTrainDefList = samplesTrainDefList,
                                    nModels = nModels,
                                    ntree = ntree)

  # Store the settings of the classifier run within the resulting object
  metaDataRun <- base::data.frame(nModels = nModels,
                            classColumn = classColumn,
                            higherClassColumn = higherClassColumn,
                            domainColumn = domainColumn,
                            maxSamplesPerType = maxSamplesPerType,
                            nANOVAgenes = nANOVAgenes,
                            nFeatures = nFeatures,
                            whichSeed = whichSeed)

  createdModelsMinority <- base::list(modelList = modelList,
                                riboModelList = riboModelList,
                                reducedFeatures = reducedFeatures,
                                metaDataRef = metaDataRef,
                                metaDataRun = metaDataRun,
                                countDataRef = countDataOG
                                )
  if (useUpsimpler == T) {
    createdModelsMinority$metaDataSynths <- synths$synthMetadataDF
    createdModelsMinority$countsSynths <- countsSynths

    relevantParametersUpsampling <- upsamplerArgs$upsimple

    #createdModelsMinority$relevantParametersUpsampling <- relevantParametersUpsampling
    if (upsamplingType == "upsimpler") {
      upsimplerParameters <- base::data.frame(nNeighbors = relevantParametersUpsampling$n_neighbors,
                                              nSamples = relevantParametersUpsampling$n_samples,
                                              excludeSamplesFromSameClass = relevantParametersUpsampling$exclude_same_class,
                                              allowDistancesBetweenDifferentClasses = relevantParametersUpsampling$allow_interclass_dists,
                                              subsampleFraction = relevantParametersUpsampling$subsample_frac,
                                              sampleWithReplacement = relevantParametersUpsampling$with_replacement,
                                              seed = relevantParametersUpsampling$rseed,
                                              scalingStrategy = relevantParametersUpsampling$scaling_strategy,
                                              scalingMaxDistance = relevantParametersUpsampling$nnmax_factor,
                                              repelFromOtherDatapoint = relevantParametersUpsampling$nnrepel_factor,
                                              negativeValues = relevantParametersUpsampling$neg_strategy
      )
    } else if (upsamplingType == "random") {
      upsimplerParameters <- data.frame(Type = "random")
    }

    createdModelsMinority$upsimplerParameters <- upsimplerParameters
    createdModelsMinority$metaDataRun$upsampleBelow <- upsampleBelow
    createdModelsMinority$metaDataRun$maxSamplesAfterUpsampling <- maxSamplesAfterUpsampling
  }

  if (saveModel == T) {
  filename <- base::paste0(directory, "/createdModelsMinority.rds")
  base::saveRDS(createdModelsMinority, file = filename)
  base::print(base::paste0("Please find the generated R-object with the created minority classification models within ", filename))
  }
  return(createdModelsMinority)
}
