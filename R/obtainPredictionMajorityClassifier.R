#' Obtain predictions from generated majority models
#'
#' This function looks at the accompanying prediction for the test samples for the different generated Majority classifier models.
#'
#' @param rotationsAndScalingsList List containing both the derived PCA-transformation information ($prList)
#' and the information needed to scale new sample input data to center the features around 0 ($scaleFeaturesList).
#' @param dataTest New samples for which the classifier should provide a classification label.
#' @param metaDataRef Metadata file containing the links between the samples and the tumor subtype diagnosis.
#' @param testSamples test sample IDs matching with samples from the test data.
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#' @param nModels How many models should be created for the Majority classifier?
#' @param vectorK Vector generated within the function 'trainKNNClassifier', calculating the optimal value for _k_ given the dataset and knn-model.
#' @return Dataframe containing the classifications for the _nModels_ different generated models, with the different folds
#' in the columns and the different samples to be predicted in the rows.
#' @import kknn
#'
obtainPredictionMajorityClassifier <- function(rotationsAndScalingsList,
                                     dataTest,
                                     metaDataRef,
                                     testSamples,
                                     classColumn,
                                     maxNeighbors = 25,
                                     nModels = 100,
                                     vectorK
) {
  for (i in base::seq(1:nModels)) {
    base::print(base::paste0("Working on model ", i))
    rotations <- rotationsAndScalingsList[["prList"]][[i]]$rotation
    varFeatures <- rotationsAndScalingsList[["scaleFeaturesList"]][[i]]$varFeatures
    meanGenes <- rotationsAndScalingsList[["scaleFeaturesList"]][[i]]$meanGenes
    sdGenes <- rotationsAndScalingsList[["scaleFeaturesList"]][[i]]$sdGenes

    # Rotate the training data
    rotatedTrainData <- rotationsAndScalingsList[["prList"]][[i]]$x %>% base::as.data.frame()
    testSamples <- base::colnames(dataTest)
    rotatedTestSamples <- base::t((dataTest[varFeatures,]-meanGenes[varFeatures])/sdGenes[varFeatures]) %*% rotations %>%
      base::as.data.frame()

    rotatedTrainData$class <- base::as.factor(metaDataRef[base::rownames(rotatedTrainData),classColumn])

    model <- kknn::kknn(class~., rotatedTrainData,
                  rotatedTestSamples,
                  distance = 1,
                  kernel = "optimal",
                  scale= T,
                  k= vectorK[i])

    prediction <- base::as.character(model$fitted.values)

    if (i == 1) {
      result <- base::data.frame(fold1 = prediction)
      base::rownames(result) <- testSamples
    } else {
      result[, base::paste0("fold", i)] <- prediction
    }
  }

  return(result)
}
