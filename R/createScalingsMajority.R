#' Create transformation information Majority Classifier
#'
#' This function is used to obtain the PCA-rotations and data scalings centered around zero
#' that are generate within the Majority Classifier. These rotations and scalings are needed
#' to transform new samples as well.
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param metaDataRef  Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Name of column in the metadata file that contains the tumor type labels.
#' @param domainColumn Name of column in the metadata file that contains the tumor domain labels.
#' @param sampleColumn Name of column in the metadata file that contains the sample identifiers.
#' @param nModels How many models should be created for the Majority classifier?
#' @param maxSamplesPerType How many samples should we maximally use per tumor (sub)type?
#' @param nComps How many principal components will be selected after PCA?
#' @param nFeatures How many of the most variable RNA-transcripts within the dataset should we select for principal component analysis (PCA)?
#' @param maxNeighbors What is the maximum number of neigbors to be used for the weighted _k_-nearest neighbor algorithm?
#' @param upsimplerModule The upsimpler Python package as imported by \code{reticulate}. If not specified, no upsampling will be performed.
#' @param upsimplerArgs List of arguments to be passed to the upsimpler algorithm, separately for the initialization of the
#' \code{Upsimpler} class ($init), as well as for the \code{upsimple} method ($upsimple).
#' @param whichSeed For reproducibility, the seed can be specified with this parameter.
#' @param outputDir Directory in which you would like to store the R-object containing the results. Default is today's date.
#' @param proteinCodingGenes What are the names of the RNA-transcripts that stand for protein-coding genes within our dataset?
#' Please supply it as a vector. This is needed for ribo-depletion correction model.
#' @param saveModel Do you want to save your generated results in an R-object (saveModel = TRUE)? Default is TRUE.
#' @return R-object containing the rotations and scalings
#' for each training data subset($rotationsAndScalingList),
#' the model to correct for the ribodepletion efficacy ($riboModelList),
#' which samples were present in each training data subset ($samplesTrainDefList),
#' which RNA-transcripts were considered for transformation in the analysis ($nonZeroGenes),
#' the metadata file associated to the reference cohort ($metaDataRef),
#' the metadata for the performed run ($metaDataRun),
#' and the matrix containing the RNA-transcript per million data that was used to create the scaling
#' and is needed for eventual classifications of new samples ($countDataRef).
#'
#' @export
#' @import caret
#'
createScalingsMajority <-  function(countDataRef,
                                    metaDataRef,
                                    classColumn,
                                    higherClassColumn,
                                    domainColumn,
                                    sampleColumn,
                                    nModels = 100,
                                    maxSamplesPerType = 50,
                                    nComps = 100,
                                    nFeatures = 2500,
                                    maxNeighbors = 25,
                                    upsimplerModule = NULL,
                                    upsimplerArgs = NULL,
                                    whichSeed = 1,
                                    outputDir = paste0("./", format(as.Date(Sys.Date(), "%Y-%m-%d"), "%Y_%m_%d")),
                                    proteinCodingGenes,
                                    saveModel = T

) {
  `%notin%` <- base::Negate(`%in%`)
  if (sampleColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the sample IDs is not present within metaDataRef. Please check the sampleColumn.")
  } else if (classColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the tumor subtype labels is not present within metaDataRef. Please check the classColumn")
  } else if (higherClassColumn %notin% base::colnames(metaDataRef)){
    base::stop("The column you specified for the tumor type labels is not present within metaDataRef. Please check the higherClassColumn")
  } else if (domainColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the tumor domain labels is not present within metaDataRef. Please check the domainColumn")
  }

  base::rownames(metaDataRef) <- metaDataRef[, sampleColumn]
  # Make sure the metadata and count data are in the right format and same order
  if (base::nrow(metaDataRef) != base::ncol(countDataRef)) {
    base::stop("The number of samples do not match between the metadata and the count data. Please make sure you include all same samples in both objects.")
  } else if (base::all(base::rownames(metaDataRef) %notin% base::colnames(countDataRef))) {
    base::stop("Your input data is not as required. Please make sure your sample IDs are stored in the sampleColumn, and in the column names of the count data")
  }


  if (base::is.numeric(countDataRef) != T) {
    base::stop("Your input data is not as required. Please make sure your countDataRef object only contains numerical count data and is a matrix.")

  }

  # Include a statement to store the classColumn, higherClassColumn and domainColumn
  base::print(base::paste0("The column used for tumor subtypes labels within the metadata, used for model training purposes, is: ", classColumn, ', containing values such as: '))
  base::print(base::unique(metaDataRef[,classColumn])[1:3])

  base::print(base::paste0("The column used for tumor type labels within the metadata, is: ", higherClassColumn,', containing values such as: '))
  base::print(base::unique(metaDataRef[,higherClassColumn])[1:3])

  base::print(base::paste0("The column used for tumor domain labels within the metadata, is: ", domainColumn, ', containing values such as: '))
  base::print(base::unique(metaDataRef[,domainColumn])[1:3])
  base::print("If any of these are incorrect, specify a different 'classColumn' (subtype), 'higherClassColumn' (tumor type) or 'domainColumn' (domain) to function as labels.")

  tumorEntitiesWithTooFewSamples <- base::table(metaDataRef[,classColumn])[base::table(metaDataRef[,classColumn]) < 3] %>% base::names()
  if (base::length(tumorEntitiesWithTooFewSamples) > 0) {
    # there are tumor entities (classes) with less than 3 samples
    # if upsimpler and its args are provided, we can upsample these classes
    # otherwise, we need to remove them
    msg <- "You have labels within your dataset that have less than 3 available samples."
    if (!is.null(upsimplerModule) && !is.null(upsimplerArgs)) {
      msg <- base::paste(msg, "Please note that these classes will be upsampled (upsimpler was provided).")
      # define the classes and the respective samples than need to be dealt with
      classesToUpsample <- tumorEntitiesWithTooFewSamples
      # select a single sample from each of these classes and make it a vector
      seedSamples <- base::lapply(classesToUpsample, function(x) {
        base::sample(metaDataRef[metaDataRef[, classColumn] == x, sampleColumn], 1)
      }) %>% base::unlist()
    } else {
      metaDataRef %<>% dplyr::filter(!!dplyr::sym(classColumn) %notin% tumorEntitiesWithTooFewSamples)
      msg <- base::paste(msg, "Please note samples with these labels have been removed.")
    }
    base::print(msg)
  }

  if (!base::dir.exists(outputDir)) {
    checkDirectory <- base::tryCatch(base::dir.create(outputDir))
    if (checkDirectory == F) {
      base::stop(base::paste0("The directory you want the classification to be saved in cannot be created due to an error in the directory path.",
                              " Please check the spelling of your specified outputDir - it is probable the parent-directory does not exist."))
    }
  }

  # Make sure you have CPM counts
  countDataRef <- base::apply(countDataRef,2,function(x) (x/base::sum(x))*1E6)

  directory <- outputDir

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
  base::print("Finished the ribocorrection")
  countDataRef <- riboModelList$counts

  # Log-transform data
  dataLogRef <- base::log(countDataRef +1)

  # Remove the genes with zero variance across the dataset

  #zeroVar <- nearZeroVar(dataCV)
  dataLogZeroVar <- base::t(dataLogRef) %>% base::as.data.frame(.)
  zeroVar <- caret::nearZeroVar(dataLogZeroVar)

  dataLogNonZero <- dataLogRef[-zeroVar, ]
  nonZeroGenes <- base::rownames(dataLogNonZero)
  # Set seed for reproducibility
  base::set.seed(whichSeed)

  # Select biomaterial IDs as training data per model

  samplesTrainDefList <- obtainTrainData(metaDataRef = metaDataRef,
                                         classColumn = classColumn,
                                         nModels = nModels,
                                         maxSamplesPerType = maxSamplesPerType)

  base::print("samplesTrainDefList created")

  samplesPerModel <- selectFeaturesMajority(dataTrain = dataLogNonZero,
                                            samplesTrainDefList = samplesTrainDefList,
                                            nFeatures = nFeatures)

  ##### UPSIMPLER INTEGRATION PART 1 START #####
  if (!is.null(upsimplerModule) && !is.null(upsimplerArgs) && length(seedSamples) > 0) {
    # things are more complicated here than in Minority: each model has a different
    # feature set. According to the upsimplerArgs$maj_strategy, we will either:
    # - "hom": homogeneous upsampling (i.e., on the union of all features for all models)
    # - "het-": heterogeneous upsampling (i.e., on the features of each model separately, before PCA)
    # - "het+": heterogeneous upsampling (i.e., on the features of each model separately, after PCA)
    if (upsimplerArgs$maj_strategy == "hom") {

      base::print("Homogeneous upsampling strategy selected")

      allMajorityFeatures <- base::lapply(samplesPerModel, base::rownames) %>% base::unlist() %>% base::unique()

      # Apply the upsimpler algorithm
      upsimplerArgs$init$dataDF <- dataLogNonZero[allMajorityFeatures,] %>% base::t() %>% base::as.data.frame()
      upsimplerArgs$init$metadataDF <- metaDataRef
      upsimplerArgs$init$class_col <- classColumn

      upsimpler <- do.call(upsimplerModule$Upsimpler, upsimplerArgs$init)

      synths <- applyUpsimpler(upsimpler = upsimpler,
                               upsimplerArgs = upsimplerArgs$upsimple,
                               targetSampleIDs = seedSamples,
                               metadataDF = metaDataRef,
                               classColumn = classColumn,
                               sampleColumn = sampleColumn)

      # distribute the synths to the models
      synthsPerModel <- obtainTrainData(metaDataRef = synths$synthMetadataDF,
                                        classColumn = classColumn,
                                        nModels = nModels,
                                        maxSamplesPerType = maxSamplesPerType)

      # map synthsPerModel from sample IDs to the actual data
      synths$synthDataDF %<>% base::t() %>% base::as.data.frame()
      for (i in base::seq_along(synthsPerModel)) {
        synthsPerModel[[i]] <- synths$synthDataDF[base::rownames(samplesPerModel[[i]]), synthsPerModel[[i]]]
      }
      synths$synthDataDF %<>% base::t() %>% base::as.data.frame()

      # add the synths to the training data
      samplesPerModel <- base::lapply(base::seq_along(samplesPerModel), function(i) {
        base::cbind(samplesPerModel[[i]], synthsPerModel[[i]])
      })

      # synths metadata is a single column; we need to add the other columns
      # in order to rbind it to the original metadata
      synths$synthMetadataDF[base::setdiff(base::names(metaDataRef), base::names(synths$synthMetadataDF))] <- NA
      metaDataRef <- base::rbind(metaDataRef, synths$synthMetadataDF)

    } else if (upsimplerArgs$maj_strategy == "het-") {

      base::print("Heterogeneous upsampling strategy selected (before feature extraction via PCA)")
      # Apply the upsimpler algorithm
      synthsPerModel <- base::list()
      for (i in base::seq_along(samplesPerModel)) {
        upsimplerArgs$init$dataDF <- samplesPerModel[[i]] %>% base::t() %>% base::as.data.frame()
        upsimplerArgs$init$metadataDF <- metaDataRef
        upsimplerArgs$init$class_col <- classColumn

        upsimpler <- do.call(upsimplerModule$Upsimpler, upsimplerArgs$init)

        synths <- applyUpsimpler(upsimpler = upsimpler,
                                 upsimplerArgs = upsimplerArgs$upsimple,
                                 targetSampleIDs = seedSamples,
                                 metadataDF = metaDataRef,
                                 classColumn = classColumn,
                                 sampleColumn = sampleColumn)

        # add the synths to the training data
        samplesPerModel[[i]] <- base::cbind(samplesPerModel[[i]], synths$synthDataDF %>% base::t() %>% base::as.data.frame())
        # append the synths metadata to the original metadata
        synths$synthMetadataDF[base::setdiff(base::names(metaDataRef), base::names(synths$synthMetadataDF))] <- NA
        # append them row by row, only if the synths are not already present in the metadata
        for (j in base::seq_len(base::nrow(synths$synthMetadataDF))) {
          if (!(base::rownames(synths$synthMetadataDF)[j] %in% base::rownames(metaDataRef))) {
            metaDataRef <- base::rbind(metaDataRef, synths$synthMetadataDF[j,])
          }
        }

        synthsPerModel[[i]] <- synths$synthDataDF
      }
    }

  }
  ####### UPSIMPLER INTEGRATION PART 1 END #######

  rotationsAndScalingsList <- getPrincipalComponents(samplesPerModel = samplesPerModel,
                                                     nComps = nComps)

  ###### UPSIMPLER INTEGRATION PART 2 START ######
  if (!is.null(upsimplerModule) && !is.null(upsimplerArgs) && length(seedSamples) > 0 && upsimplerArgs$maj_strategy == "het+") {
    base::print("Heterogeneous upsampling strategy selected (after feature extraction via PCA)")
    # Apply the upsimpler algorithm
    synthsPerModel <- base::list()
    for (i in base::seq_along(rotationsAndScalingsList$prList)) {
      upsimplerArgs$init$dataDF <- rotationsAndScalingsList$prList[[i]]$x %>% base::as.data.frame()
      upsimplerArgs$init$metadataDF <- metaDataRef
      upsimplerArgs$init$class_col <- classColumn

      # in PCA space, negative values are not a problem
      upsimplerArgs$upsimple$neg_strategy <- NULL

      upsimpler <- do.call(upsimplerModule$Upsimpler, upsimplerArgs$init)

      synths <- applyUpsimpler(upsimpler = upsimpler,
                               upsimplerArgs = upsimplerArgs$upsimple,
                               targetSampleIDs = seedSamples,
                               metadataDF = metaDataRef,
                               classColumn = classColumn,
                               sampleColumn = sampleColumn)

      # add the synths to the training data
      rotationsAndScalingsList$prList[[i]]$x <- base::rbind(rotationsAndScalingsList$prList[[i]]$x, synths$synthDataDF)
      # append the synths metadata to the original metadata
      synths$synthMetadataDF[base::setdiff(base::names(metaDataRef), base::names(synths$synthMetadataDF))] <- NA
      # append them row by row, only if the synths are not already present in the metadata
      for (j in base::seq_len(base::nrow(synths$synthMetadataDF))) {
        if (!(base::rownames(synths$synthMetadataDF)[j] %in% base::rownames(metaDataRef))) {
          metaDataRef <- base::rbind(metaDataRef, synths$synthMetadataDF[j,])
        }
      }

      synthsPerModel[[i]] <- synths$synthDataDF
    }
  }
  ####### UPSIMPLER INTEGRATION PART 2 END #######

  # Store the settings of the classifier run within the resulting object
  metaDataRun <- base::data.frame(nModels = nModels,
                            classColumn = classColumn,
                            higherClassColumn = higherClassColumn,
                            domainColumn = domainColumn,
                            sampleColumn = sampleColumn,
                            maxSamplesPerType = maxSamplesPerType,
                            maxNeighbors = maxNeighbors,
                            nFeatures = nFeatures,
                            nComps = nComps,
                            whichSeed = whichSeed)

  createdModelsMajority <- base::list(rotationsAndScalingsList = rotationsAndScalingsList,
                                riboModelList = riboModelList,
                                nonZeroGenes = nonZeroGenes,
                                metaDataRef = metaDataRef,
                                metaDataRun = metaDataRun)

  # if upsampling was performed, add the synths to the result
  if (!is.null(upsimplerModule) && !is.null(upsimplerArgs) && length(seedSamples) > 0) {
    createdModelsMajority$synthsPerModel <- synthsPerModel
  }

  if (saveModel == T) {

    filename <- base::paste0(directory, "/createdModelsMajority.rds")

    base::saveRDS(createdModelsMajority, file = filename)
    base::print(base::paste0("Please find the generated R-object with the created majority models and PCA-transformations within ", filename))
  }

  return(createdModelsMajority)
}
