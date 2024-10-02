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
                                    whichSeed = 1,
                                    outputDir = paste0("./", format(as.Date(Sys.Date(), "%Y-%m-%d"), "%Y_%m_%d")),
                                    proteinCodingGenes,
                                    saveModel = T

) {
  #countDataOG <- countDataRef
  `%notin%` <- base::Negate(`%in%`)

  checkFormatInputData(sampleColumn = sampleColumn,
                       classColumn = classColumn,
                       higherClassColumn = higherClassColumn,
                       domainColumn = domainColumn,
                       metaDataRef = metaDataRef,
                       countDataRef = countDataRef,
                       outputDir = outputDir
                       )
  base::rownames(metaDataRef) <- metaDataRef[, sampleColumn]

  tumorEntitiesWithTooFewSamples <- base::table(metaDataRef[,classColumn])[base::table(metaDataRef[,classColumn]) < 3] %>% base::names()
  if (base::length(tumorEntitiesWithTooFewSamples) >0) {

    metaDataRef %<>% dplyr::filter(!!dplyr::sym(classColumn) %notin% tumorEntitiesWithTooFewSamples)
    base::cat("You have labels within your dataset that have less than 3 available samples.  Please note samples with these labels have been removed.")
    #stop("You have tumor subtypes within your dataset that have less than 3 available samples. Please remove all tumor types with too few samples. ")

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
  base::cat("\nFinished the ribocorrection\n")
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

  base::cat("\nsamplesTrainDefList created\n")

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

  # Store the settings of the classifier run within the resulting object
  metaDataRun <- base::data.frame(nModels = nModels,
                            classColumn = classColumn,
                            higherClassColumn = higherClassColumn,
                            domainColumn = domainColumn,
                            maxSamplesPerType = maxSamplesPerType,
                            maxNeighbors = maxNeighbors,
                            nFeatures = nFeatures,
                            nComps = nComps,
                            whichSeed = whichSeed)

  createdModelsMajority <- base::list(rotationsAndScalingsList = rotationsAndScalingsList,
                                riboModelList = riboModelList,
                                samplesTrainDefList = samplesTrainDefList,
                                nonZeroGenes = nonZeroGenes,
                                metaDataRef = metaDataRef,
                                metaDataRun = metaDataRun,
                                countDataRef = countDataOG
                                )

  if (saveModel == T) {

    filename <- base::paste0(directory, "/createdModelsMajority.rds")

    base::saveRDS(createdModelsMajority, file = filename)
    base::cat(base::paste0("\nPlease find the generated R-object with the created majority models and PCA-transformations within ", filename, "\n"))
  }

  return(createdModelsMajority)
}
