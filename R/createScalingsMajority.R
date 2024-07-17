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
  countDataOG <- countDataRef
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
  if (base::length(tumorEntitiesWithTooFewSamples) >0) {

    metaDataRef %<>% dplyr::filter(!!dplyr::sym(classColumn) %notin% tumorEntitiesWithTooFewSamples)
    base::print("You have labels within your dataset that have less than 3 available samples.  Please note samples with these labels have been removed.")
    #stop("You have tumor subtypes within your dataset that have less than 3 available samples. Please remove all tumor types with too few samples. ")

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


  rotationsAndScalingsList <- getPrincipalComponents(dataTrain = dataLogNonZero,
                                                     samplesTrainDefList,
                                                     nModels = nModels,
                                                     nFeatures = nFeatures,
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
                                metaDataRun = metaDataRun#,
                                #countDataRef = countDataOG
                                )

  if (saveModel == T) {

    filename <- base::paste0(directory, "/createdModelsMajority.rds")

    base::saveRDS(createdModelsMajority, file = filename)
    base::print(base::paste0("Please find the generated R-object with the created majority models and PCA-transformations within ", filename))
  }

  return(createdModelsMajority)
}
