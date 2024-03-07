#' Get PCA transformation Majority Classifier
#'
#' Function to find the optimal principal component analysis (PCA) transformation of the reference cohort
#' so that a maximal variance within the data is explained. With the parameters, the number of variable features
#' and the eventual total number of principal components can be adjusted.
#'
#' @param dataTrain Data to be used for training, containing the whole reference cohort.
#' @param samplesTrainDefList List of different training data subsets used for the majority voting system.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param nModels How many models should be created for the majority voting system?
#' @param nFeatures How many of the most variable genes within the dataset should we select for principal component analysis (PCA)?
#' @param nComps How many principal components will be selected after PCA?
#'
#' @return List containing both the derived PCA-transformation information ($prList)
#' and the information needed to scale new sample input data to center the features around 0 ($scaleFeaturesList).
#'
getPrincipalComponents <- function(dataTrain,
         samplesTrainDefList,
         classColumn,
         nModels = 100,
         nFeatures = 2500,
         nComps = 100
) {

  # The genes are in the columns right now, so the variance should be determined per column (removing)

  prList <- list()
  scaleFeaturesList <- list()
  for (i in seq(1:nModels)) {
    print(paste("Working on model", i))
    samplesTrainDef <- samplesTrainDefList[[i]]

    train.data <- dataTrain[, samplesTrainDef]

    varGenes <- apply(train.data,1 ,var)

    varFeatures <- names(varGenes)[order(varGenes,decreasing = T)][c(1:nFeatures)]

    dataTrainFiltered <- as.data.frame(train.data) %>%
      filter(rownames(.) %in% varFeatures)

    meanGenes <- apply(dataTrainFiltered, 1, mean,
                       na.rm = T)


    sdGenes <- apply(dataTrainFiltered, 1, sd)

    dataScale <- apply(dataTrainFiltered[varFeatures,], 2, function(x) (x-meanGenes[varFeatures])/sdGenes[varFeatures])

    pr <- stats::prcomp(t(dataScale),
                 rank.=nComps)

    prList[[i]] <- pr
    scaleFeatures <- list(varFeatures = varFeatures,
                          meanGenes = meanGenes,
                          sdGenes = sdGenes)
    scaleFeaturesList[[i]] <- scaleFeatures

    rotationsAndScalingsList <- list(prList = prList,
                                     scaleFeaturesList = scaleFeaturesList)
  }
  return(rotationsAndScalingsList)
}
