#' Get PCA transformation Majority Classifier
#'
#' Function to find the optimal principal component analysis (PCA) transformation of the reference cohort
#' so that a maximal variance within the data is explained. With the parameters, the eventual total number
#' of principal components can be adjusted.
#'
#' @param samplesPerModel List of different training data sample subsets used for subsetting the available training dataset.
#' These contain the actual samples, not just the IDs (as they were selected by \code{selectFeaturesMajority} or
#' as returned by \code{upsimplerResampling}).
#' @param nComps How many principal components will be selected after PCA?
#'
#' @return List containing both the derived PCA-transformation information ($prList)
#' and the information needed to scale new sample input data to center the features around 0 ($scaleFeaturesList).
#'
getPrincipalComponents <- function(samplesPerModel, nComps = 100) {
  prList <- base::list()
  scaleFeaturesList <- base::list()
  for (i in base::seq_along(samplesPerModel)) {
    base::print(base::paste("Working on model", i))

    train.data <- samplesPerModel[[i]]
    varFeatures <- base::rownames(train.data)

    meanGenes <- base::apply(train.data, 1, base::mean,
                       na.rm = T)

    sdGenes <- base::apply(train.data, 1, stats::sd)

    dataScale <- base::apply(train.data, 2, function(x) (x - meanGenes) / sdGenes)

    pr <- stats::prcomp(base::t(dataScale), rank.=nComps)

    prList[[i]] <- pr
    scaleFeatures <- base::list(varFeatures = varFeatures,
                                meanGenes = meanGenes,
                                sdGenes = sdGenes)

    scaleFeaturesList[[i]] <- scaleFeatures

    rotationsAndScalingsList <- base::list(prList = prList,
                                           scaleFeaturesList = scaleFeaturesList)
  }
  return(rotationsAndScalingsList)
}
