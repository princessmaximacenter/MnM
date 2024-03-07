#' Obtain predictions from generated majority models
#'
#' This function looks at the accompanying prediction for the test samples for the different generated Majority Classifier models.
#'
#' @param rotationsAndScalingsList List containing both the derived PCA-transformation information ($prList)
#' and the information needed to scale new sample input data to center the features around 0 ($scaleFeaturesList).
#' @param dataTrain Data to be used for training, containing the whole reference cohort.
#' @param dataTest Data for which the classifier should predict a label.
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param samplesTrainDefList List of different training data subsets used for the majority voting system.
#' @param testSamples IDs for test samples from the test data.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param maxNeighbours What is the maximum number of neigbours to be used for the weighted _k_-nearest neighbor algorithm?
#' @param nModels How many models should be created for the majority voting system?
#'
#' @return Dataframe containing the predictions for the _nModels_ different generated models, with the different folds
#' in the columns and the different samples to be predicted in the rows.
#' @export
#' @import magrittr
#'
obtainPredictionMajorityClassifier <- function(rotationsAndScalingsList,
                                     dataTrain,
                                     dataTest,
                                     metaDataRef,
                                     samplesTrainDefList,
                                     testSamples,
                                     classColumn = classColumn,
                                     maxNeighbours = 25,
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
                         scale=F,ks=c(1:maxNeighbours))

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
