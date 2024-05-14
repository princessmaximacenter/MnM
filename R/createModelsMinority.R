#' Generate Minority classifier models
#'
#' Setup to obtain the models for the Minority classifier.
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param metaDataRef Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param domainColumn Column in the metadata file that contains the tumor domain labels.
#' @param sampleColumn Column in the metadata file that contains the sample identifiers.
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
#'
#'
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
  } else if (base::all(base::rownames(metaDataRef) %notin% colnames(countDataRef))) {
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

  countDataRef <- countDataRef[, base::rownames(metaDataRef)]
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
  countDataRef <- riboModelList$counts

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
                                metaDataRun = metaDataRun)

  if (saveModel == T) {
  filename <- base::paste0(directory, "/createdModelsMinority.rds")
  base::saveRDS(createdModelsMinority, file = filename)
  base::print(base::paste0("Please find the generated R-object with the created minority classification models within ", filename))
  }
  return(createdModelsMinority)
}
