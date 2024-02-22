combineSeedPredictions <- function(nSeeds,
         minorityDir,
         majorityDir,
         subtype = F,
         metaDataRef,
         higherClassColumn,
         classColumn) {
for (i in seq(1:nSeeds)) {
  minorityDoc <- paste0(minorityDir, "seed",i, "/crossValidationMinorityResults.rds")
  majorityDoc <- paste0(majorityDir, "seed",i, "/crossValidationMajorityResults.rds")
  minority <- readRDS(minorityDoc)
  majority <- readRDS(majorityDoc)

  probabilitiesMinority <- obtainProbabilities(minority,
                                               crossValidation = T,
                                               nModels = 100
  )
  probabilitiesMajority <- obtainProbabilities(majority,
                                               crossValidation = T,
                                               nModels = 100)

  if (subtype == F) {
  linkClassAndHigherClass <- minority$metaData[ , c(classColumn, higherClassColumn)] %>% unique

  probabilitiesMinority <- changeSubtypeNameToType(probabilitiesMinority,
                                                   linkClassAndHigherClass = linkClassAndHigherClass,
                                                   classColumn = classColumn,
                                                   higherClassColumn = higherClassColumn)
  probabilitiesMajority <- changeSubtypeNameToType(probabilitiesMajority,
                                                   linkClassAndHigherClass = linkClassAndHigherClass,
                                                   classColumn = classColumn,
                                                   higherClassColumn = higherClassColumn)
  }
  MMProbabilityList <- getMMProbabilities(minorityProbability = probabilitiesMinority,
                                          majorityProbability = probabilitiesMajority)

  keys <- names(MMProbabilityList)

  if (i == 1) {
    MMProbabilityListFinal <- MMProbabilityList

  } else {
    MMProbabilityListFinal <- setNames(mapply(c, MMProbabilityListFinal[keys], MMProbabilityList[keys]), keys)
  }
}
MMProbabilityListFinalFinal <- MMProbabilityListFinal
for (i in seq(1:length(MMProbabilityListFinal))) {
  MMProbabilityListFinalFinal[[i]] <- tapply(MMProbabilityListFinal[[i]], names(MMProbabilityListFinal[[i]]), sum)
}

MMProbabilityListFinalFinal <- lapply(MMProbabilityListFinalFinal, function(x) x / 10)

if (subtype == T) {
predictionsMMFinal <- getMajorityPredictions(minority = minority,
                                             MMProbabilityList = MMProbabilityListFinalFinal,
                                             higherClassColumn = classColumn,
                                             crossValidation = T,
                                             metaDataRef = metaDataRef,
                                             subtype = T)
} else {
  predictionsMMFinal <- getMajorityPredictions(minority = minority,
                                               MMProbabilityList = MMProbabilityListFinalFinal,
                                               higherClassColumn = higherClassColumn,
                                               crossValidation = T,
                                               metaDataRef = metaDataRef,
                                               subtype = F)

}
predictionsList <- list(predictionsMMFinal = predictionsMMFinal,
                        MMProbabilityList = MMProbabilityListFinalFinal)

return(predictionsList)

}
