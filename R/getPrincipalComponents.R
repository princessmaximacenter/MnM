#' Get PCA transformation Majority Classifier
#'
#' Function to find the optimal principal component analysis (PCA) transformation of the reference cohort
#' so that a maximal variance within the data is explained. With the parameters, the number of variable features
#' and the eventual total number of principal components can be adjusted.
#'
#' @param dataTrain Data to be used for model training, containing the all samples available within the training dataset.
#' @param samplesTrainDefList List of different training data sample subsets used for subsetting the available training dataset (dataTrain).
#' @param nModels How many models should be created for the classifier?
#' @param nFeatures How many of the most variable RNA-transcripts within the dataset should we select for principal component analysis (PCA)?
#' @param nComps How many principal components will be selected after PCA?
#'
#' @return List containing both the derived PCA-transformation information ($prList)
#' and the information needed to scale new sample input data to center the features around 0 ($scaleFeaturesList).
#'
getPrincipalComponents <- function(dataTrain,
                                   scaleFeaturesList,
         samplesTrainDefList,
         nModels = 100,
         nComps = 100
) {

  # The RNA-transcripts are in the columns right now, so the variance should be determined per column

  prList <- base::list()
  for (i in base::seq(1:nModels)) {
    base::cat(base::paste("Working on model", i, "\n"))
    samplesTrainDef <- samplesTrainDefList[[i]]

    meanGenes <- scaleFeaturesList[[i]]$meanGenes

    sdGenes <- scaleFeaturesList[[i]]$sdGenes

    varFeatures <- scaleFeaturesList[[i]]$varFeatures

    train.data <- dataTrain[, samplesTrainDef]

    dataTrainFiltered <- base::as.data.frame(train.data) %>%
      dplyr::filter(base::rownames(.) %in% varFeatures)

    dataScale <- base::apply(dataTrainFiltered[varFeatures,], 2, function(x) (x-meanGenes[varFeatures])/sdGenes[varFeatures])

    pr <- stats::prcomp(base::t(dataScale),
                 rank.=nComps)

    prList[[i]] <- pr

    rotationsAndScalingsList <- base::list(prList = prList,
                                     scaleFeaturesList = scaleFeaturesList)
  }
  return(rotationsAndScalingsList)
}
