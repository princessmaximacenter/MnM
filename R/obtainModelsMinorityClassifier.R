#' Title
#'
#' @param dataTrain
#' @param samplesTrainDefList
#' @param nModels
#' @param ntree
#'
#' @return
#' @export
#'
obtainModelsMinorityClassifier <- function(dataTrain,
                                           samplesTrainDefList,
                                           nModels,
                                           ntree = 500
) {
  modelList <- list()

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
