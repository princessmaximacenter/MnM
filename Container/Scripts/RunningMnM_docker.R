###########
#Dedicated R Script to running the MnM Classifier via containerization
#Date 06-12-2024
#Project: MnM Classifier
###########

"Loading Packages"
library("devtools")
library("yaml")
library("optparse") # Library for parsing command-line options
#Load MnM Package
library("MnM")
library("remotes")

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              help = "YAML input filename")
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the input filename is provided
if (is.null(opt$input)) {
  stop("Error: Input YAML filename is missing. Use -i or --input flag to specify the filename.")
}

"Load YAML File for Variable Use (Via Command Line)"
MnM_config <- read_yaml(opt$input)

"Loading Metadata"
##Metadata for Reference Cohort
metaDataRef <- read.table(MnM_config$Inputs$reference_metadata, sep = "\t", row.names = 1)

##Metadata for Test Cohort
metaDataTest <- read.table(MnM_config$Inputs$test_metadata, sep = "\t", row.names = 1)

"Loading Count Data"
##Reference Cohort
countDataRef <- read.table(MnM_config$Inputs$reference_count_data, sep = "\t",
                           row.names = 1) %>% as.matrix()

##Test Cohort
countDataTest <- read.table(MnM_config$Inputs$test_count_data, sep = "\t",
                         row.names = 1) %>% as.matrix()

"Data Prep"
## Remove Samples that have less than 3 available samples (if applicable)
filter_samples = MnM_config$modelsVariables$filter_samples

if (filter_samples){
  metaDataFiltered <- metaDataRef %>% group_by(tumorSubtype) %>%
    filter(n() >= 3) %>% ungroup()

  countDataFiltered <- countDataRef[,metaDataFiltered$sampleID]

  ## Check whether the same number of samples are present within the metadata and the count data.
  ncol(countDataFiltered) == nrow(metaDataFiltered)


} else {
  #Reassign Column Values
  countDataRef <- countDataRef[, metaDataRef$sampleID]
  ncol(countDataRef) == nrow(metaDataRef)
}

"#Performing classification for new samples"
#Input Protein Coding Genes
proteinCodingGenes <- read.csv(MnM_config$modelsVariables$proteinCodingGenes, header = T) %>%
  as.vector() %>% deframe

#Set Variable for model training to true or false
model_training = MnM_config$modelsVariables$model_training


#update values according to R yaml file
if (model_training) {
  ##createModelsMinority dedicated to generate Minority predictions
  modelsMinority <- createModelsMinority(
    countDataRef = countDataRef,
    metaDataRef = metaDataRef,
    classColumn = MnM_config$modelsMinorityVariables$classColumn,
    higherClassColumn = MnM_config$modelsMinorityVariables$higherClassColumn,
    domainColumn = MnM_config$modelsMinorityVariables$domainColumn,
    sampleColumn = MnM_config$modelsMinorityVariables$sampleColumn,
    nModels = MnM_config$modelsTrainingVariables$nModels,
    outputDir = MnM_config$modelsTrainingVariables$outputDir,
    proteinCodingGenes = proteinCodingGenes
  )

  ##createScalingsMajority dedicated to generate Majority predictions
  modelsMajority <- createScalingsMajority(
    countDataRef = countDataRef,
    metaDataRef = metaDataRef,
    sampleColumn = MnM_config$modelsMajorityVariables$sampleColumn,
    classColumn = MnM_config$modelsMajorityVariables$classColumn,
    higherClassColumn = MnM_config$modelsMajorityVariables$higherClassColumn,
    domainColumn = MnM_config$modelsMajorityVariables$domainColumn,
    nModels = MnM_config$modelsVariables$nModels,
    outputDir = MnM_config$modelsVariables$outputDir,
    proteinCodingGenes = proteinCodingGenes
  )

} else {
  # Import Models if training is not necessary
  modelsMinority <- readRDS(MnM_config$inputs$minority_model)
  modelsMajority <- readRDS(MnM_config$inputs$majority_model)
}

##Revise Predictions per Minority and Majority Classifiers
#Add conditional when training model again on new reference cohort
predictionsTestMinority <- newPredictionsMinority(
                      createdModelsMinority = modelsMinority,
                      countDataNew = countDataTest,
                      outputDir = MnM_config$modelsVariables$outputDir)

predictionsTestMajority <- newPredictionsMajority(
                      createdModelsMajority = modelsMajority,
                      countDataNew = countDataTest,
                      outputDir = MnM_config$modelsVariables$outputDir,
                      countDataRef = countDataRef)

"Predictions with Integrated Minority and Majority Classifiers"
##Integrate the Minority and Majority classifications, set CV to F
predictionsMMTestList <- integrateMM(minority = predictionsTestMinority,
                        majority = predictionsTestMajority,
                        subtype = T)

"Obtain Final Predictions"
##Call predictions from IntegrateMM
predictionsMMTest <- predictionsMMTestList$predictionsMMFinal

##Match diagnostic label to samples
predictionsMMTest$originalCall <- metaDataTest[rownames(predictionsMMTest),"tumorSubtype"]

"Export Final Predictions Table"
# Add biomaterial_id column using row names
predictionsMMTest$biomaterial_id <- rownames(predictionsMMTest)
# Reorder columns
predictionsMMTest <- predictionsMMTest %>%
  select(biomaterial_id, everything())

write.table(predictionsMMTest, file.path(MnM_config$modelsVariables$outputDir, "predictionsMMTest.tsv"),
            sep = "\t", row.names = F)

"Export Predictions and Models as R Objects to the Outputs Directory"
#Save Predictions
saveRDS(predictionsTestMinority, file = file.path(MnM_config$modelsVariables$outputDir, "predictionsTestMinority.rds"))
saveRDS(predictionsTestMajority, file = file.path(MnM_config$modelsVariables$outputDir, "predictionsTestMajority.rds"))
saveRDS(predictionsMMTestList, file = file.path(MnM_config$modelsVariables$outputDir, "predictionsMMTestList.rds"))


#Save Trained Models
if(model_training) {
  saveRDS(modelsMajority, file = file.path(MnM_config$modelsVariables$outputDir, "modelsMinority.rds"))
  saveRDS(modelsMajority, file = file.path(MnM_config$modelsVariables$outputDir, "modelsMajority.rds"))
}