#' Generate Models Minority Classifier
#'
#' Function to generated weighted Random Forest (RF) models to predict tumor (sub)type labels for new patient samples.
#' Each model is generated using a different subset of the training data, as specified in the _samplesTrainDefList_.
#'
#' @param dataTrain Data to be used for training, containing the whole reference cohort.
#' @param samplesTrainDefList List of different training data subsets used for the majority voting system.
#' @param nModels How many models should be created for the majority voting system?
#' @param ntree How many trees should we use during the weighted RF procedure?
#'
#' @return List containing the different Minority Classifier models to use for subsequent tumor (sub)type predictions.
#'
obtainModelsMinorityClassifier <- function(dataTrain,
                                           samplesTrainDefList,
                                           nModels,
                                           ntree = 500
) {
  modelList <- list()

  #classesVal <- as.character(dataTrain$class) %>% table()
  #probabilityClasses <- 1 / classesVal
  #classwt <- as.numeric(probabilityClasses)
  for (i in seq(1:nModels)) {
    print(paste("Working on model", i))
    trainSamples <- samplesTrainDefList[[i]]

    train.data <- dataTrain[rownames(dataTrain) %in% trainSamples,]

    train.category <- as.character(train.data$class)
    train.data %<>% select(-c("class"))

    # Determine class weights for the weighted RF
    classesVal <- table(train.category)
    probabilityClasses <- 1/classesVal
    classwt <- as.numeric(probabilityClasses)

    # Generate RF model on training subset
    model <- randomForest(x = train.data,
                          y = as.factor(train.category),
                          importance = T,
                          ntree = ntree,
                          proximity = F,
                          classwt = classwt)

    modelList[[i]] <- model
  }

  return(modelList)
}
