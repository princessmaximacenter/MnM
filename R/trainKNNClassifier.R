#' Calculate optimal value for k in a given model
#'
#' @param rotationsAndScalingsList List containing both the derived PCA-transformation information ($prList)
#' and the information needed to scale new sample input data to center the features around 0 ($scaleFeaturesList).
#' @param metaDataRef Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#' @param maxNeighbors What is the maximum number of neighbors to be used for the weighted _k_-nearest neighbor algorithm?
#' @param nModels How many models should be created for the Majority classifier?
#' @import kknn
#'
#' @return Vector with optimized values for k within the different weighted kNN models
#'
trainKNNClassifier <- function(rotationsAndScalingsList,
                               metaDataRef,
                               classColumn,
                               maxNeighbors = 25,
                               nModels = 100) {

vectorK <- c()
for (i in base::seq(1:nModels)) {
  base::cat(base::paste0("Working on model ", i, "\n"))

  # Rotate the training data
  rotatedTrainData <- rotationsAndScalingsList[["prList"]][[i]]$x %>% base::as.data.frame()
  # Add training class labels
  rotatedTrainData$class <- base::as.factor(metaDataRef[base::rownames(rotatedTrainData),classColumn])

  rotatedTrainDataK <- rotatedTrainData[base::grep("\\.",base::rownames(rotatedTrainData),invert=T),]

  # Calculate optimal value for k for the given data in the specific feature space
  if (nrow(rotatedTrainDataK) < 30) {
    nfolds <- ceiling(nrow(rotatedTrainDataK) / 10)
  } else {
    nfolds <- 10
  }

  if (maxNeighbors > nrow(rotatedTrainDataK)) {
    maxNeighborsNew <- nrow(rotatedTrainDataK)
  } else {
    maxNeighborsNew <- maxNeighbors
  }

  kTrain <- kknn::train.kknn(class~., rotatedTrainDataK,
                             distance = 1,
                             kernel = "optimal",
                             scale=F,ks=c(1:maxNeighborsNew),
                             kcv= nfolds)

  # Store calculated optimal k-value
  vectorK[i] <- kTrain$best.parameters$k
}

return(vectorK)
}
