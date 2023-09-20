calculateMulticlassROC <- function(predictionsMM,
                                   minority,
                                   majority,
                                   crossValidation,
                                   nModels
                                   ) {

  probabilityForAUC <- data.frame(
    originalCall = predictionsMM$originalCall)

  rownames(probabilityForAUC) <- rownames(predictionsMM)

  tumors <- unique(predictionsMM$originalCall)

  for (i in seq(1:length(tumors))) {
    probabilityForAUC[, tumors[i]] <- 0

  }

  probabilitiesMinority <- obtainProbabilities(minority,
                                               crossValidation = crossValidation,
                                               nModels = nModels
  )
  probabilitiesMajority <- obtainProbabilities(majority,
                                               crossValidation = crossValidation,
                                               nModels = nModels)


  if (changeNames == T) {
    for (i in seq(1:length(probabilitiesMinority))) {
      for (j in seq(1:length(substituteNames))) {
        names(probabilitiesMinority[[i]]) <- gsub(substituteNames[j], substituteBy, names(probabilitiesMinority[[i]]))

      }
      probabilitiesMinority[[i]] <- tapply(probabilitiesMinority[[i]], names(probabilitiesMinority[[i]]), sum)
    }

    for (i in seq(1:length(probabilitiesMajority))) {
      for (j in seq(1:length(substituteNames))) {
        names(probabilitiesMajority[[i]]) <- gsub(substituteNames[j], substituteBy, names(probabilitiesMajority[[i]]))

      }
      probabilitiesMajority[[i]] <- tapply(probabilitiesMajority[[i]], names(probabilitiesMajority[[i]]), sum)
    }

  }

  for (i in seq(1:length(MMProbabilityList))) {
    predictions <- MMProbabilityList[[i]]
    patient <- names(MMProbabilityList[i])
    for (j in seq(1:length(predictions))) {
      probabilityForAUC[patient, names(predictions)[j]] <- as.numeric(predictions[j])
    }
  }



}
