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
#'
#' @return R-object containing the generated RF-models ($modelList), the model for the ribodepletion correction ($riboModelList),
#'  the features that were eventually used for the weighted RF within the different folds ($reducedFeaturesList),
#'   the metadata file associated to the reference cohort ($metaData)
#'  and the metadata for the performed run ($metaDataRun).
#' @export
#'
createModelsMinority <-  function(countDataRef,
                                  metaDataRef,
                                  meanExpression = 5,
                                  classColumn,
                                  nModels = 100,
                                  nANOVAgenes = 1000,
                                  maxSamplesPerType = NA,
                                  ntree = 500,
                                  howManyFeatures = 300,
                                  whichSeed = 1,
                                  outputDir = "~/Documents/data/output/",
                                  proteinDir = "~/surfdrive/Shared/Kemmeren group/Research_Projects/RNA_classification_FW/data/input/"

) {

  # Make sure you have CPM counts
  countDataRef <- apply(countDataRef,2,function(x) (x/sum(x))*1E6)


  directory <- paste0(outputDir, format(as.Date(Sys.Date(), "%Y-%m-%d"), "%m_%d_%Y"), "/")
  if (!dir.exists(directory)) {
    dir.create(directory) }

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
                                    howManyFeatures = howManyFeatures)

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

  filename <- paste0(directory, "/createdModelsMinority.rds")
  write_rds(createdModelsMinority, file = filename)
  return(createdModelsMinority)
}
