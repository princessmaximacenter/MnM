#' Combine the results of different cross-validation runs
#'
#' To check for algorithm stability, cross-validation setups can be run multiple times
#' to see whether results remain similar. This function combines the results of different runs into one,
#' coming with a final prediction for each sample used during the cross-validation setup.
#'
#' @param nSeeds How many seeds was the cross-validation setup run with?
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param subtype Do you want to combine the classifications on the subtype level?
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param nModels How many models were created for the majority voting system?
#'
#' @return List containing the a dataframe with the final predictions (predictionsMMFinal),
#' and a list with the probability scores for all classification labels that were assigned to samples (MMProbabilityList).
#'
#'
#' predictionsMMFinal contains the top 3 final predictions ($predict{2,3}) with their accompanying probability scores
#' ($probability{1,2,3} and the original diagnosis label ($originalCall).
#'
#'
combineSeedPredictions <- function(nSeeds,
                                   nModels,
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
                                               nModels = nModels
  )
  probabilitiesMajority <- obtainProbabilities(majority,
                                               crossValidation = T,
                                               nModels = nModels)

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

MMProbabilityListFinalFinal <- lapply(MMProbabilityListFinalFinal, function(x) x / nSeeds)

if (subtype == T) {
predictionsMMFinal <- getMajorityPredictions(minority = minority,
                                             MMProbabilityList = MMProbabilityListFinalFinal,
                                             higherClassColumn = classColumn,
                                             crossValidation = T,
                                             metaDataRef = metaDataRef,
                                             subtype = subtype)
} else {
  predictionsMMFinal <- getMajorityPredictions(minority = minority,
                                               MMProbabilityList = MMProbabilityListFinalFinal,
                                               higherClassColumn = higherClassColumn,
                                               crossValidation = T,
                                               metaDataRef = metaDataRef,
                                               subtype = subtype)

}
predictionsList <- list(predictionsMMFinal = predictionsMMFinal,
                        MMProbabilityList = MMProbabilityListFinalFinal)

return(predictionsList)

}
