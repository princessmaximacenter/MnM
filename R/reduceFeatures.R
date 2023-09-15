#' Reduce features Minority Classifier with RF-importance
#'
#' This functions' purpose is to further reduce the number of features within the Minority Classifier.
#' Ten Random Forest (RF) models will be generated using the reference cohort training data,
#' from which the importance scores for the different genes will be extracted.
#' The _howManyFeatures_ genes with the highest importance scores will be selected for the final Minority Classifier training.
#'
#' @param dataTrain Data to be used for training, containing the whole reference cohort.
#' @param samplesTrainDefList List of different training data subsets used for the majority voting system.
#' @param ntree How many trees should we use during the weighted Random Forest (RF) procedure?
#' @param nModels How many models should be created for the majority voting system?
#' @param howManyFeatures How many features should we keep after determining
#' the most important genes using the Random Forest Importance Score?
#'
#' @return Names of the top most important genes present in the training data,
#' containing most information to distinguish the different tumor (sub)types.
#'
reduceFeatures <- function(dataTrain,
                           samplesTrainDefList,
                           ntree = 500,
                           nModels =10,
                           howManyFeatures = 300) {

  modelList <- list()
  nModels <- min(nModels, 10)
  for (i in seq(1:nModels)) {
    samplesTrainDef <- samplesTrainDefList[[i]]

    train.data <- dataTrain[rownames(dataTrain) %in% samplesTrainDef,]

    train.category <- as.character(train.data$class)

    train.data %<>% select(-c("class"))

    classesVal <- table(train.category)
    probabilityClasses <- 1/classesVal

    classwt <- as.numeric(probabilityClasses)

    model <- randomForest(x = train.data, y = as.factor(train.category),
                          importance = T, ntree = ntree,
                          proximity = F, classwt = classwt)
    modelList[[i]] <- model
  }

  accuracyValuesDF <- createFeatureDF(modelList=modelList,
                                      whichAccuracyMeasure = "MeanDecreaseAccuracy",
                                      nANOVAgenes=nANOVAgenes)

  meanAccuracyValuesDF <- apply(accuracyValuesDF, 1, mean)

  topFeatures <- meanAccuracyValuesDF %>%
    sort(decreasing = T) %>%
    head(n = howManyFeatures)

  topFeaturesNamesAccuracy <- names(topFeatures)

  return(topFeaturesNamesAccuracy)
}
