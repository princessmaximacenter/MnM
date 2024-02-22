#' Create transformation information Majority Classifier
#'
#' This function is used to obtain the PCA-rotations and data scalings around zero
#' that are generate within the Majority Classifier. These rotations and scalings are needed
#' to transform new samples as well.
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data. Patients are in the columns,
#' different genes in the rows.
#' @param metaDataRef Metadata file containing the links between the patients and
#' the tumor (sub)type diagnosis within the reference cohort.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param nModels How many models should be created for the majority voting system?
#' @param maxSamplesPerType How many samples should we maximally use per tumor (sub)type?
#' @param nComps How many principal components will be selected after PCA?
#' @param nFeatures How many of the most variable genes within the dataset should we select for principal component analysis (PCA)?
#' @param maxNeighbours What is the maximum number of neigbours to be used for the weighted _k_-nearest neighbor algorithm?
#' @param whichSeed For reproducibility, the seed can be specified with this parameter.
#' @param outputDir Directory in which you would like to store the R-object containing the results.
#' @param proteinDir In which directory can we find the file specifying the names of protein-coding genes within our dataset?
#' @param patientColumn Column in the metadata file that contains the patient labels.
#' @param saveModel
#' @return R-object containing the rotations and scalings
#' for each reference cohort subset($rotationsAndScalingList),
#' the model to correct for the ribodepletion efficacy ($riboModelList),
#' which samples were present in each subset ($samplesTrainDefList),
#' which genes were considered for transformation in the analysis ($nonZeroGenes),
#' the metadata file associated to the reference cohort ($metaData),
#' and metadata for the performed run ($metaDataRun).
#'
#' @export
#'
createScalingsMajority <-  function(countDataRef,
                                    metaDataRef,
                                    classColumn,
                                    patientColumn,
                                    nModels = 100,
                                    maxSamplesPerType = 50,
                                    nComps = 100,
                                    nFeatures = 2500,
                                    maxNeighbours = 25,
                                    whichSeed = 1,
                                    outputDir = "./",
                                    proteinDir,
                                    saveModel = T,
                                    saveRiboModels = F

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
                                       outputDir = directory,
                                       saveRiboModels = saveRiboModels
    )

  } else {
    riboModelList <- readRDS(riboCountFile)
  }
  print("Finished the ribocorrection")
  countDataRef <- riboModelList$counts

  # Log-transform data
  dataLogRef <- log(countDataRef +1)

  # Remove the genes with zero variance across the dataset

  #zeroVar <- nearZeroVar(dataCV)
  dataLogZeroVar <- t(dataLogRef) %>% as.data.frame(.)
  zeroVar <- nearZeroVar(dataLogZeroVar)

  dataLogNonZero <- dataLogRef[-zeroVar, ]
  nonZeroGenes <- rownames(dataLogNonZero)
  # Set seed for reproducibility
  set.seed(whichSeed)

  # Select biomaterial IDs as training data per model

  samplesTrainDefList <- obtainTrainData(metaDataRef = metaDataRef,
                                         classColumn = classColumn,
                                         nModels = nModels,
                                         maxSamplesPerType = maxSamplesPerType)

  print("samplesTrainDefList created")


  rotationsAndScalingsList <- getPrincipalComponents(dataTrain = dataLogNonZero,
                                                     samplesTrainDefList,
                                                     classColumn = classColumn,
                                                     nModels = nModels,
                                                     nFeatures = nFeatures,
                                                     nComps = nComps
  )

  # Store the settings of the classifier run within the resulting object
  metaDataRun <- data.frame(nModels = nModels,
                            classColumn = classColumn,
                            maxSamplesPerType = maxSamplesPerType,
                            nFeatures = nFeatures,
                            nComps = nComps,
                            whichSeed = whichSeed)

  createdModelsMajority <- list(rotationsAndScalingsList = rotationsAndScalingsList,
                                riboModelList = riboModelList,
                                samplesTrainDefList = samplesTrainDefList,
                                nonZeroGenes = nonZeroGenes,
                                metaData = metaDataRef,
                                metaDataRun = metaDataRun)

  if (saveModel == T) {
    filename <- paste0(directory, "/createdModelsMajority.rds")
    write_rds(createdModelsMajority, file = filename)
  }

  return(createdModelsMajority)
}
