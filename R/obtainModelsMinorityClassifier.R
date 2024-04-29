#' Generate models Minority classifier
#'
#' Function to generated weighted Random Forest (RF) models to predict tumor (sub)type labels for new patient samples.
#' Each model is generated using a different subset of the training data, as specified in the _samplesTrainDefList_.
#'
#' @param dataTrain Data to be used for model training, containing the all samples available within the training dataset.
#' @param samplesTrainDefList List of different training data sample subsets used for subsetting the available training dataset (dataTrain).
#' @param nModels How many models should be created for the classifier?
#' @param ntree How many trees should we generate during the weighted RF procedure?
#'
#' @return List containing the different Minority classifier models to use for subsequent tumor (sub)type predictions.
#' @import randomForest
obtainModelsMinorityClassifier <- function(dataTrain,
                                           samplesTrainDefList,
                                           nModels = 100,
                                           ntree = 500
) {
  modelList <- list()

  for (i in seq(1:nModels)) {
    print(paste("Working on model", i))
    trainSamples <- samplesTrainDefList[[i]]

    train.data <- dataTrain[rownames(dataTrain) %in% trainSamples,]

    train.category <- as.character(train.data$class)
    train.data %<>% dplyr::select(-c("class"))

    # Determine class weights for the weighted RF
    classesVal <- base::table(train.category)
    probabilityClasses <- 1/classesVal
    classwt <- as.numeric(probabilityClasses)

    # Generate RF model on training subset
    model <- randomForest::randomForest(x = train.data,
                          y = as.factor(train.category),
                          importance = T,
                          ntree = ntree,
                          proximity = F,
                          classwt = classwt)

    modelList[[i]] <- model
  }

  return(modelList)
}
