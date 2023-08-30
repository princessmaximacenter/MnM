#' Title
#'
#' @param rotationsAndScalingsList
#' @param dataTrain data to be used for training, containing the whole dataset.
#' @param dataTest data for which the classifier should predict a label.
#' @param metaData information on the labels of the training data.
#' @param samplesTrainDefList list of different training data subsets used for the majority voting system.
#' @param testSamples
#' @param classColumn
#' @param nComps
#' @param maxNeighbours
#' @param nModels
#'
#' @return
#' @export
#'
#' @examples
obtainPredictionMajority <- function(rotationsAndScalingsList,
                                     dataTrain,
                                     dataTest,
                                     metaData,
                                     samplesTrainDefList,
                                     testSamples,
                                     classColumn = classColumn,
                                     nComps = 100,
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

    rotatedTestSamples <- t((dataTest[varFeatures,]-meanGenes[varFeatures])/sdGenes[varFeatures]) %*% rotations %>%
      as.data.frame()

    rotatedTrainData$class <- as.factor(metaData[rownames(rotatedTrainData),classColumn])

    rotatedTrainDataK <- rotatedTrainData[grep("\\.",rownames(rotatedTrainData),invert=T),]

    kTrain <- train.kknn(class~., rotatedTrainDataK,
                         distance = 1,
                         kernel = "optimal",
                         scale=F,ks=c(1:maxNeighbours))

    model <- kknn(class~., rotatedTrainData,
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
