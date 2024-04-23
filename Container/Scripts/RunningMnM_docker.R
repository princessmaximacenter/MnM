###########
#Dedicated R Script to running the MnM Classifier via docker/singularity
#Date 23-04-2024
#Project Tumor-Classification
###########

"Loading Packages"
library("devtools")
library("yaml")
library("optparse") # Library for parsing command-line options
#Load MnM Package
libray("MnM")

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

# Load YAML file for variable use
MnM_config <- read_yaml(opt$input)

"Load YAML File for Variable Use"
MnM_config <- read_yaml("/app/Inputs/RunningMnM_variables_docker_inputs.yaml")

"Loading Metadata"
##Metadata for Reference Cohort
metaDataRef <- read.table(MnM_config$inputs$reference_metadata, sep = "\t", row.names = 1)

##Metadata for Test Cohort
metaDataTest <- read.table(MnM_config$inputs$test_metadata, sep = "\t", row.names = 1)

"Loading Count Data"
##Reference Cohort
countDataRef <- read.table(MnM_config$inputs$reference_count_data, sep = "\t",
                           row.names = 1) %>% as.matrix()

##Test Cohort
testCounts <- read.table(MnM_config$inputs$test_count_data, sep = "\t",
                         row.names = 1) %>% as.matrix()

"Data Prep"
##Assign biomaterial_id values
rownames(metaDataRef) <- metaDataRef$biomaterial_id

## Align matching rows between countDataRef and metaDataRef
countDataRef <- countDataRef[, rownames(metaDataRef)]

"Temporary Filtering on B-ALL Data to test model (will be removed later)"
##Select the samples within the metadata from tumor type B-ALL
selectedSamples <- metaDataRef %>% filter(Disease_sub_class == "B-ALL") %>% rownames()

##Generate a filtered metadata and count data object with only the selected samples.
metaDataFiltered <- metaDataRef[selectedSamples,] %>% as.data.frame()

countDataFiltered <- countDataRef[,selectedSamples]

## Check whether the same number of samples are present within the metadata and the count data.
ncol(countDataFiltered) == nrow(metaDataFiltered)

"#Performing classification for new samples"
model_training = MnM_config$model_training_variables$model_training


#update values according to R ya
if (model_training) {
  ##createModelsMinority dedicated to generate Minority predictions
  modelsMinority <- createModelsMinority(
    countDataRef = countDataFiltered,
    metaDataRef = metaDataFiltered,
    patientColumn = "biomaterial_id",
    meanExpression = 5,
    classColumn = "Disease_sub_specification1",
    nModels = 10,
    nANOVAgenes = 1000,
    maxSamplesPerType = 3,
    ntree = 500,
    howManyFeatures = 300,
    whichSeed = 1,
    outputDir = "./Outputs/",
    proteinDir = "./Inputs/",
    saveModel = T
  )

  ##createScalingsMajority dedicated to generate Majority predictions
  modelsMajority <- createScalingsMajority(
    countDataRef = countDataFiltered,
    metaDataRef = metaDataFiltered,
    patientColumn = "biomaterial_id",
    classColumn = "Disease_sub_specification1",
    nModels = 10,
    maxSamplesPerType = 50,
    nComps = 100,
    nFeatures = 2500,
    maxNeighbours = 25,
    whichSeed = 1,
    outputDir = "./Outputs/",
    proteinDir =  "./Inputs/",
    saveModel = T
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
                      countDataNew = testCounts,
                      outputDir = "./Outputs/",
                      saveModel = T)

predictionsTestMajority <- newPredictionsMajority(
                      createdModelsMajority = modelsMajority,
                      countDataNew = testCounts,
                      outputDir = "./Outputs/",
                      saveModel = T,
                      countDataRef = countDataFiltered,
                      classColumn = "Disease_sub_specification1",
                      nModels = 10) #Default is 100

"Predictions with Integrated Minority and Majority Classifiers"
##Integrate the Minority and Majority classifications, set CV to F
predictionsMMTestList <- integrateMM(minority = predictionsTestMinority,
                        predictionsTestMajority,
                        metaDataRef = metaDataFiltered,
                        nModels = 10, #Default 100
                        subtype = T,
                        classColumn = "Disease_sub_specification1",
                        higherClassColumn = "Disease_sub_class",
                        crossValidation = F)

"Obtain Final Predictions"
##Call predictions from IntegrateMM
predictionsMMTest <- predictionsMMTestList$predictionsMMFinal

##Match diagnostic label to samples
predictionsMMTest$originalCall <- metaDataTest[rownames(predictionsMMTest),"Disease_sub_specification1"]

"Export Final Predictions Table"
# Add biomaterial_id column using row names
predictionsMMTest$biomaterial_id <- rownames(predictionsMMTest)
# Reorder columns
predictionsMMTest <- predictionsMMTest %>%
  select(biomaterial_id, everything())

write.table(predictionsMMTest, "./Outputs/predictionsMMTest.tsv", sep = "\t", row.names = F)

"Export Predictions and Models as R Objects to the Outputs Directory"
#Save Predictions
saveRDS(predictionsTestMinority, file = "./Outputs/predictionsTestMinority.rds")
saveRDS(predictionsTestMajority, file = "./Outputs/redictionsTestMajority.rds")
saveRDS(predictionsMMTestList, file = "./Outputs/predictionsMMTestList.rds")

#Save Trained Models
if(model_training) {
  saveRDS(modelsMinority, file = "./Outputs/modelsMinority.rds")
  saveRDS(modelsMajority, file = "./Outputs/modelsMajority.rds")
}
