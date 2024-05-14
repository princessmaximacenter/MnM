#' Reduce features Minority classifier with Random Forest importance score
#'
#' This functions' purpose is to further reduce the number of features within the Minority classifier.
#' Ten Random Forest (RF) models will be generated using the reference cohort training data,
#' from which the importance scores for the different RNA-transcripts will be extracted.
#' Each model is generated using a different subset of the training data, as specified in the _samplesTrainDefList_.
#' The _nFeatures_ RNA-transcripts with the highest importance scores will be selected for the final Minority Classifier training.
#'
#' @param dataTrain Data to be used for model training, containing the all samples available within the training dataset.
#' @param samplesTrainDefList List of different training data sample subsets used for subsetting the available training dataset (dataTrain).
#' @param ntree How many trees should we generate during the weighted RF procedure?
#' @param nModels How many models should be created for the majority voting system?
#' @param nFeatures How many features should we keep after determining
#' the most important RNA-transcripts using the RF Importance Score?
#' @param nANOVAgenes How many RNA-transcripts did we select during the ANOVA procedure?
#' @import randomForest
#' @return Names of the top most important RNA-transcripts present in the training data,
#' containing most information to distinguish the different tumor (sub)types.

reduceFeatures <- function(dataTrain,
                           samplesTrainDefList,
                           ntree = 500,
                           nModels =10,
                           nFeatures = 300,
                           nANOVAgenes) {

  modelList <- list()
  nModels <- base::min(nModels, 100)
  for (i in base::seq(1:nModels)) {
    samplesTrainDef <- samplesTrainDefList[[i]]

    train.data <- dataTrain[base::rownames(dataTrain) %in% samplesTrainDef,]

    train.category <- base::as.character(train.data$class)

    train.data %<>% dplyr::select(-c("class"))

    classesVal <- base::table(train.category)
    probabilityClasses <- 1/classesVal

    classwt <- base::as.numeric(probabilityClasses)

    model <- randomForest::randomForest(x = train.data, y = as.factor(train.category),
                          importance = T, ntree = ntree,
                          proximity = F, classwt = classwt)
    modelList[[i]] <- model
  }

  accuracyValuesDF <- createFeatureDF(modelList=modelList,
                                      whichAccuracyMeasure = "MeanDecreaseAccuracy",
                                      nANOVAgenes=nANOVAgenes)

  meanAccuracyValuesDF <- base::apply(accuracyValuesDF, 1, base::mean)

  topFeatures <- meanAccuracyValuesDF %>%
    base::sort(decreasing = T) %>%
    utils::head(n = nFeatures)

  topFeaturesNamesAccuracy <- base::names(topFeatures)

  return(topFeaturesNamesAccuracy)
}
