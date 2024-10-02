#' Extract most variable features from training dataset
#'
#' @param dataTrain Data to be used for model training, containing the all samples available within the training dataset.
#' @param samplesTrainDefList List of different training data sample subsets used for subsetting the available training dataset (dataTrain).
#' @param nFeatures How many of the most variable RNA-transcripts within the dataset should we select for principal component analysis (PCA)?
#' @param nModels How many models should be created for the classifier?
#'
#' @return List containing the most variable features, portraying information needed to scale new sample input data to center the features around 0.
#'
getVarFeaturesMajority <- function(
    dataTrain,
    samplesTrainDefList,
    nFeatures = 2500,
    nModels = 100

    ) {


  scaleFeaturesList <- base::list()
  for (i in base::seq(1:nModels)) {
    base::cat(base::paste("Working on model", i, "\n"))
    samplesTrainDef <- samplesTrainDefList[[i]]

    train.data <- dataTrain[, samplesTrainDef]

    varGenes <- base::apply(train.data,1 ,stats::var)

    varFeatures <- base::names(varGenes)[base::order(varGenes,decreasing = T)][c(1:nFeatures)]

    dataTrainFiltered <- base::as.data.frame(train.data) %>%
      dplyr::filter(base::rownames(.) %in% varFeatures)

    meanGenes <- base::apply(dataTrainFiltered, 1, base::mean,
                             na.rm = T)

    sdGenes <- base::apply(dataTrainFiltered, 1, stats::sd)

    scaleFeatures <- base::list(varFeatures = varFeatures,
                                meanGenes = meanGenes,
                                sdGenes = sdGenes)

    scaleFeaturesList[[i]] <- scaleFeatures

  }

  return(scaleFeaturesList)

}
