
# MnM Classifier via Containerization

## Overview

The **RunningMnM.r** script is dedicated to executing the MnM (Minority and Majority) Classifier for tumor classification. It integrates metadata and count data for reference and test cohorts, performs classification using pre-trained or newly trained models, and outputs predictions and saved models.

This script uses containerization for streamlined execution, ensuring reproducibility and consistency in results. The configuration inputs are provided via a **YAML** file.

## Key Features

- **Command-line Input**: Accepts a YAML configuration file as input to dynamically specify metadata, count data, and model training parameters.
- **Data Preprocessing**: Filters samples based on criteria and aligns metadata with count data.
- **Model Training and Prediction**:
  - **Minority Classifier**: Trains models on minority subtypes.
  - **Majority Classifier**: Trains scaling models for majority subtypes.
- **Integration**: Combines predictions from Minority and Majority classifiers.
- **Export Options**: Saves predictions and trained models to specified output directories.

---

## Contents

### Script: RunningMnM.r

#### 1. **Package Loading**
   Loads necessary R packages such as `optparse`, `yaml`, and `MnM` for command-line parsing, configuration management, and classification.

#### 2. **Input Parsing**
   Uses `optparse` to read the YAML configuration file specified by the `-i` or `--input` flag. Throws an error if the input is not provided.

#### 3. **Configuration Management**
   Reads the YAML file to extract input paths, model variables, and other runtime parameters.

#### 4. **Data Loading**
   Imports metadata and count data for the reference and test cohorts using paths specified in the YAML file.

#### 5. **Data Preparation**
   - Filters reference samples if the `filter_samples` option is enabled.
   - Ensures alignment between metadata and count data.

#### 6. **Model Training and Predictions**
   - If `model_training` is enabled:
     - Trains Minority and Majority classifiers using `createModelsMinority` and `createScalingsMajority` functions.
   - Otherwise:
     - Loads pre-trained models from the specified paths.
   - Generates predictions for test cohorts using `newPredictionsMinority` and `newPredictionsMajority`.
   - Integrates results using the `integrateMM` function.

#### 7. **Output Generation**
   - Exports final predictions as a TSV file.
   - Saves trained models and prediction objects in RDS format.

---

## Configuration File: YAML Input

### Example Input Structure
```yaml
# Script Configuration
script:
  name: "MnM Classifier via Docker/Singularity"
  date: "11-12-2024"
  project: "Tumor-Classification"

# File Paths
Inputs:
  reference_metadata: "/app/Inputs/ClassifierMetadataReferenceBALL.csv"
  test_metadata: "/app/Inputs/ClassifierMetadataTestBALL.csv"
  reference_count_data: "/app/Inputs/countDataReferenceBALL.csv"
  test_count_data: "/app/Inputs/countDataTestBALL.csv"
  majority_model: "/app/SavedModels/modelsMajority.rds"
  minority_model: "/app/SavedModels/modelsMinority.rds"

# Model Training values
modelsVariables:
  model_training: TRUE
  filter_samples: FALSE
  whichSeed: 1
  outputDir: "/app/Outputs"
  proteinCodingGenes: "/app/Inputs/proteinCodingGenesGencode31.csv"
  saveModel: TRUE
  saveRiboModels: FALSE
  nModels: 10

# Specific values for createModelsMinority
modelsMinorityVariables:
  sampleColumn: "sampleID"
  classColumn: "tumorSubtype"
  higherClassColumn: "tumorType"
  domainColumn: "Domain"
  nANOVAgenes: 1000
  meanExpression: 5
  maxSamplesPerType: 3
  ntree: 500
  howManyFeatures: 300

# Specific values for createScalingsMajority
modelsMajorityVariables:
  sampleColumn: "sampleID"
  classColumn: "tumorSubtype"
  higherClassColumn: "tumorType"
  domainColumn: "Domain"
  maxSamplesPerType: 50
  nComps: 100
  nFeatures: 2500
  maxNeighbours: 25
  whichSeed: 1

# Specific values for Predictions
classifierPredictions:
  outputDir: "/app/Outputs"
  saveModel: "TRUE"
```

### Key Sections
1. **File Paths**:
   - Specify paths for metadata, count data, and pre-trained models.
2. **Model Variables**:
   - Control training parameters such as the number of models, filtering options, and output directories.
3. **Minority and Majority Classifier Variables**:
   - Define specific parameters for training models.

---

## Building the Docker Image

To build the MnM Docker container, use the provided Dockerfile. Ensure Docker is installed and accessible on your system. Run the following command:

```bash
docker build -t mnm:1.0.0 .
```

Here:
- `mnm` is the image name.
- `1.0.0` is the version tag.

The build process installs all system dependencies, R packages, and the MnM library from GitHub.

---

## Running the MnM Classifier via Docker

To execute the MnM classifier using Docker:

```bash
docker run --rm -v $(pwd):/app mnm:1.0.0 Rscript RunningMnM.r -i /app/config.yaml
```

- `--rm`: Ensures the container is removed after execution.
- `-v $(pwd):/app`: Mounts the current directory to `/app` inside the container, enabling access to input files and scripts.
- `Rscript RunningMnM.r`: Executes the classifier script.

---

## Outputs

- **Predictions Table**:
  - TSV file containing final predictions with associated metadata.
- **Saved Models**:
  - RDS files of trained models and prediction objects for future use.

---

## Dependencies

The Docker image includes the following:

### System Dependencies
- Build tools (`build-essential`)
- Development libraries for SSL, Curl, XML, SQLite, etc.
- Graphics libraries (e.g., Cairo, Fontconfig, Freetype)

### R Version
- **R 4.4.0** (Base image: `r-base:4.4.0`)

### R Packages
- Core utilities: `yaml`, `optparse`, `devtools`
- MnM library installed from GitHub (`princessmaximacenter/MnM`, `dev` branch)
- Statistical and machine learning tools: `randomForest`, `glmnet`, `caret`, etc.

---

## Notes

- Ensure that input files are formatted correctly (e.g., tab-delimited for metadata and count data).
- Review YAML configuration for accuracy before execution.

---

## Contact
For further assistance, please contact the genomicscore@prinsesmaximacentrum.nl