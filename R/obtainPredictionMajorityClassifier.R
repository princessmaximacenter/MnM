#' Obtain predictions from generated majority models
#'
#' This function looks at the accompanying prediction for the test samples for the different generated Majority classifier models.
#'
#' @param rotationsAndScalingsList List containing both the derived PCA-transformation information ($prList)
#' and the information needed to scale new sample input data to center the features around 0 ($scaleFeaturesList).
#' @param dataTrain Data to be used for training, containing the all samples available within the training dataset.
#' @param dataTest New samples for which the classifier should provide a classification label.
#' @param metaDataRef Metadata file containing the links between the samples and the tumor subtype diagnosis.
#' @param testSamples test sample IDs matching with samples from the test data.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param maxNeighbors What is the maximum number of neighbors to be used for the weighted _k_-nearest neighbor algorithm?
#' @param nModels How many models should be created for the Majority classifier?
#'
#' @return Dataframe containing the classifications for the _nModels_ different generated models, with the different folds
#' in the columns and the different samples to be predicted in the rows.
#' @import kknn
#'
obtainPredictionMajorityClassifier <- function(rotationsAndScalingsList,
                                     dataTrain,
                                     dataTest,
                                     metaDataRef,
                                     testSamples,
                                     classColumn,
                                     maxNeighbors = 25,
                                     nModels = 100
) {
  for (i in seq(1:nModels)) {
    print(paste0("Working on model ", i))
    rotations <- rotationsAndScalingsList[["prList"]][[i]]$rotation
    varFeatures <- rotationsAndScalingsList[["scaleFeaturesList"]][[i]]$varFeatures
    meanGenes <- rotationsAndScalingsList[["scaleFeaturesList"]][[i]]$meanGenes
    sdGenes <- rotationsAndScalingsList[["scaleFeaturesList"]][[i]]$sdGenes

    # Rotate the training data
    rotatedTrainData <- rotationsAndScalingsList[["prList"]][[i]]$x %>% as.data.frame()
    testSamples <- colnames(dataTest)
    rotatedTestSamples <- t((dataTest[varFeatures,]-meanGenes[varFeatures])/sdGenes[varFeatures]) %*% rotations %>%
      as.data.frame()

    rotatedTrainData$class <- as.factor(metaDataRef[rownames(rotatedTrainData),classColumn])

    rotatedTrainDataK <- rotatedTrainData[grep("\\.",rownames(rotatedTrainData),invert=T),]

    kTrain <- kknn::train.kknn(class~., rotatedTrainDataK,
                         distance = 1,
                         kernel = "optimal",
                         scale=F,ks=c(1:maxNeighbors))

    model <- kknn::kknn(class~., rotatedTrainData,
                  rotatedTestSamples,
                  distance = 1,
                  kernel = "optimal",
                  scale= T,
                  k= kTrain$best.parameters$k)

    prediction <- as.character(model$fitted.values)

    if (i == 1) {
      result <- data.frame(fold1 = prediction)
      rownames(result) <- testSamples
    } else {
      result[, paste0("fold", i)] <- prediction
    }
  }

  return(result)
}
