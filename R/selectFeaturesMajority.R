#' Select the most variable features for each Majority Classifier model
#'
#' @param dataTrain Data to be used for model training, containing all the samples available within the training dataset.
#' @param samplesTrainDefList List of different training data sample subsets used for subsetting the available
#' training dataset (dataTrain).
#' @param nFeatures How many of the most variable RNA-transcripts within the dataset should we select for the
#' Majority Classifier? (default = 2500)
#'
#' @return List with the same length as the number of models, comprising dataframes with the selected samples
#' per model
#'
selectFeaturesMajority <- function(dataTrain,
                                   samplesTrainDefList,
                                   nFeatures = 2500) {

  selectedSamplesPerModel <- list()
  for (i in base::seq_along(samplesTrainDefList)) {
    base::print(base::paste("Feature selection for model", i))

    samplesTrainDef <- samplesTrainDefList[[i]]

    train.data <- dataTrain[, samplesTrainDef]

    varGenes <- base::apply(train.data, 1, stats::var)

    varFeatures <- base::names(varGenes)[base::order(varGenes, decreasing = T)][c(1:nFeatures)]

    dataTrainFiltered <- base::as.data.frame(train.data) %>%
      dplyr::filter(base::rownames(.) %in% varFeatures)

    selectedSamplesPerModel[[i]] <- dataTrainFiltered
  }

  return(selectedSamplesPerModel)
}
