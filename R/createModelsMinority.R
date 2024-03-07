#' Generate Minority classifier models
#'
#' Setup to obtain the models for the Minority classifier.
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data. Patients are in the columns,
#' different genes in the rows.
#' @param metaDataRef Reference cohort metadata file containing the links between the patients and the tumor (sub)type diagnosis.
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
#' @param patientColumn Column in the metadata file that contains the patient labels.
#' @param saveModel Do you want to save your generated RF-models in an R-object?
#'
#' @return R-object containing the generated RF-models ($modelList), the model for the ribodepletion correction ($riboModelList),
#'  the features that were eventually used for the weighted RF within the different folds ($reducedFeaturesList),
#'   the metadata file associated to the reference cohort ($metaData)
#'  and the metadata for the performed run ($metaDataRun).
#'
#'  @import tidyverse dplyr magrittr randomForest caret glmnet
#' @export
#'
createModelsMinority <-  function(countDataRef,
                                  metaDataRef,
                                  patientColumn,
                                  meanExpression = 5,
                                  classColumn,
                                  nModels = 100,
                                  nANOVAgenes = 1000,
                                  maxSamplesPerType = 3,
                                  ntree = 500,
                                  howManyFeatures = 300,
                                  whichSeed = 1,
                                  outputDir = "./",
                                  proteinDir,
                                  saveModel = T

) {

  `%notin%` <- Negate(`%in%`)
  rownames(metaDataRef) <- metaDataRef[, patientColumn]
  # Make sure the metadata and count data are in the right format and same order
  if (nrow(metaDataRef) != ncol(countDataRef)) {
    stop("The number of samples do not match between the metadata and the count data. Please make sure you include all same samples in both objects.")
  } else if (all(rownames(metaDataRef) %notin% colnames(countDataRef))) {
    stop("Your input data is not as required. Please make sure your patient IDs are within the row names of the metadata, and in the column names of the count data")
  }

  if (is.numeric(countDataRef) != T) {
    stop("Your input data is not as required. Please make sure your countDataRef object only contains numerical count data and is a matrix.")

  }
  tumorEntitiesWithTooFewSamples <- table(metaDataRef[,classColumn])[table(metaDataRef[,classColumn]) < 3] %>% names()
  if (length(tumorEntitiesWithTooFewSamples) >0) {

    metaDataRef %<>% filter(!!sym(classColumn) %notin% tumorEntitiesWithTooFewSamples)
    print("You have labels within your dataset that have less than 3 available samples.  Please note samples with these labels have been removed.")
    #stop("You have tumor subtypes within your dataset that have less than 3 available samples. Please remove all tumor types with too few samples. ")

  }
  countDataRef <- countDataRef[, rownames(metaDataRef)]
  # Make sure you have CPM counts
  countDataRef <- apply(countDataRef,2,function(x) (x/sum(x))*1E6)


  directory <- paste0(outputDir, format(as.Date(Sys.Date(), "%Y-%m-%d"), "%m_%d_%Y"), "/")

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

  print("We will now start with the selection of ANOVA-genes")
  # Set seed for reproducibility
  set.seed(whichSeed)


  # Run an ANOVA to select the top n genes from the training data for use in the further classification process
  interestingANOVAgenes <- selectAnovaGenes(metaDataRef = metaDataRef,
                                         countDataRef = countDataRef,
                                         nANOVAgenes = nANOVAgenes, # How many ANOVA genes
                                         classColumn = classColumn)

  print("We have finished with the ANOVA gene selection")
  # Select the ANOVA genes within the log-transformed data
  dataLogRef <- dataLogRef[ ,interestingANOVAgenes]

  # Select biomaterial IDs as training data per model
  set.seed(whichSeed)
  samplesTrainDefList <- obtainTrainData(metaDataRef = metaDataRef,
                                         classColumn = classColumn,
                                         nModels = nModels,
                                         maxSamplesPerType = maxSamplesPerType)

  # Add class label to the dataset
  dataLogRef$class <- as.character(metaDataRef[rownames(dataLogRef),classColumn])

  print("Starting to reduce features")
  # Reduce features using RF feature importance for accuracy
  set.seed(whichSeed)
  reducedFeatures <- reduceFeatures(dataTrain = dataLogRef,
                                    samplesTrainDefList = samplesTrainDefList,
                                    ntree = 500,
                                    nModels = nModels,
                                    howManyFeatures = howManyFeatures,
                                    nANOVAgenes = nANOVAgenes)

  dataLogRef <- dataLogRef[,c(reducedFeatures, "class")]

  print("Initiating RF")

  # Start the modelling of the data within the different compositions of training data
  set.seed(whichSeed)
  modelList <- obtainModelsMinorityClassifier(dataTrain = dataLogRef,
                                    samplesTrainDefList = samplesTrainDefList,
                                    nModels = nModels,
                                    ntree = ntree)

  # Store the settings of the classifier run within the resulting object
  metaDataRun <- data.frame(nModels = nModels,
                            classColumn = classColumn,
                            maxSamplesPerType = maxSamplesPerType,
                            nANOVAgenes = nANOVAgenes,
                            howManyFeatures = howManyFeatures,
                            whichSeed = whichSeed)

  createdModelsMinority <- list(modelList = modelList,
                                riboModelList = riboModelList,
                                reducedFeatures = reducedFeatures,
                                metaData = metaDataRef,
                                metaDataRun = metaDataRun)

  if (saveModel == T) {
  filename <- paste0(directory, "/createdModelsMinority.rds")
  if (!dir.exists(directory)) {
    dir.create(directory) }
  saveRDS(createdModelsMinority, file = filename)
  }
  return(createdModelsMinority)
}
