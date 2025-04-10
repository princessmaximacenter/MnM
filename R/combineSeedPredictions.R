#' Combine the results of different cross-validation runs
#'
#' To check for algorithm stability, cross-validation setups can be run multiple times
#' to see whether results remain similar. This function combines the results of different runs into one,
#' coming with a final prediction for each sample used during the cross-validation setup.
#'
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param subtype Do you want to combine the classifications on the subtype level (subtype = TRUE)?

#' @return List containing the a dataframe with the final classifications (predictionsMMFinal),
#' a list with the probability scores for all classification labels that were assigned to samples (MMProbabilityList),
#' the metadata used for the creation of the Minority and Majority models ($metaDataRef),
#' and the chosen settings of the Minority classifier ($metaDataRun)
#' to keep track of the column for the tumor subtype, tumor type and domain labels.
#'
#' predictionsMMFinal contains the top 3 final classifications ($predict{2,3}) with their accompanying probability scores
#' ($probability{1,2,3}) and the original diagnosis label ($originalCall).
#'
#' MMProbabilityList is a list containing all samples as individual list entries,
#' with their different probability scores per label.
#'
#' @export
combineSeedPredictions <- function(
         minorityDir,
         majorityDir,
         subtype = F
         ) {

  allDirsMinority <- base::list.dirs(minorityDir, recursive = F)
  allDirsMajority <- base::list.dirs(majorityDir, recursive = F)
  selectedDirsMinority <- allDirsMinority[base::grep("seed", allDirsMinority)]
  selectedDirsMajority <- allDirsMajority[base::grep("seed", allDirsMajority)]

  if (base::length(selectedDirsMinority) != base::length(selectedDirsMajority)) {
    base::stop("The number of models for the minority and majority classifier are not the same.
         Please check your models within the minorityDir and majorityDir")
  } else if (!base::identical(base::sub(majorityDir, "", selectedDirsMajority), base::sub(minorityDir, "", selectedDirsMinority)) ) {
    base::stop(base::paste("It seems that classifications from the minority and majority classifier have not been run using the same seed.",
        "\nPlease make sure you run the crossvalidation with the same seed for complementary classifications."))
  }

for (i in base::seq(1:base::length(selectedDirsMajority))) {
  minorityDoc <- base::paste0(selectedDirsMinority[i],"/crossValidationMinorityResults.rds")
  majorityDoc <- base::paste0(selectedDirsMajority[i],"/crossValidationMajorityResults.rds")
  minority <- base::readRDS(minorityDoc)
  majority <- base::readRDS(majorityDoc)

  if (i == 1) {
    classColumn <- minority$metaDataRun$classColumn
    higherClassColumn <- minority$metaDataRun$higherClassColumn
  }

  probabilitiesMinority <- obtainProbabilities(minority)
  probabilitiesMajority <- obtainProbabilities(majority)

  if (subtype == F) {
  linkClassAndHigherClass <- minority$metaDataRef[ , c(classColumn, higherClassColumn)] %>% base::unique()

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

  keys <- base::names(MMProbabilityList)

  if (i == 1) {
    MMProbabilityListFinal <- MMProbabilityList

  } else {
    MMProbabilityListFinal <- stats::setNames(base::mapply(c, MMProbabilityListFinal[keys], MMProbabilityList[keys]), keys)
  }
}
MMProbabilityListFinalFinal <- MMProbabilityListFinal
for (i in base::seq(1:base::length(MMProbabilityListFinal))) {
  MMProbabilityListFinalFinal[[i]] <- base::tapply(MMProbabilityListFinal[[i]], base::names(MMProbabilityListFinal[[i]]), base::sum)
}

MMProbabilityListFinalFinal <- base::lapply(MMProbabilityListFinalFinal, function(x) x / base::length(selectedDirsMinority))

if (subtype == T) {
predictionsMMFinal <- getTopClassifications(minority = minority,
                                             MMProbabilityList = MMProbabilityListFinalFinal,
                                             higherClassColumn = classColumn,
                                             subtype = subtype)
} else {
  predictionsMMFinal <- getTopClassifications(minority = minority,
                                               MMProbabilityList = MMProbabilityListFinalFinal,
                                               higherClassColumn = higherClassColumn,
                                               subtype = subtype)

}
predictionsList <- base::list(predictionsMMFinal = predictionsMMFinal,
                        MMProbabilityList = MMProbabilityListFinalFinal,
                        metaDataRef = minority$metaDataRef,
                        metaDataRun = minority$metaDataRun)

return(predictionsList)

}
