#' Title
#'
#' @param minority R-object that contains the results from the minority classifier
#' @param majority R-object that contains the results from the majority classifier
#' @param metaDataRef metadata-file for the reference cohort
#' @param classColumn column name within metadata-file that contains the cancer subtype-labels
#' @param higherClassColumn column name within metadata-file that contains the cancer type labels
#' @param crossValidation specify whether the results are from the cross-validation setup or not
#'
#' @return
#' @export
#'
#' @examples
integrateMM <- function(minority,
                        majority,
                        metaDataRef,
                        classColumn,
                        higherClassColumn,
                        crossValidation = T) {

  probabilitiesMinority <- obtainProbabilities(minority, crossValidation = crossValidation)
  probabilitiesMajority <- obtainProbabilities(majority, crossValidation = crossValidation)
  linkClassAndHigherClass <- metaDataRef[ , c(classColumn, higherClassColumn)] %>% unique

 # Concatenate not malignant blood samples into one 'tumor type' category

    notMalignant <- c("Not malignant bone marrow",
                      "Not malignant blood",
                      "Bone marrow failure"
    )
    gsubName <- "Not malignant (Hemato)"



  for (i in seq(1:length(notMalignant))) {
    linkClassAndHigherClass$Disease_sub_class <- gsub(notMalignant[i], gsubName, linkClassAndHigherClass$Disease_sub_class)

  }

  minorityProbabilityTumorType <- changeSubtypeNameToType(probabilitiesMinority,
                                                          linkClassAndHigherClass = linkClassAndHigherClass)
  majorityProbabilityTumorType <- changeSubtypeNameToType(probabilitiesMajority,
                                                          linkClassAndHigherClass = linkClassAndHigherClass)

  MMProbabilityList <- getMMProbabilities(minorityProbabilityTumorType = minorityProbabilityTumorType,
                                          majorityProbabilityTumorType = majorityProbabilityTumorType)

  predictionsMM <- getMajorityPredictions(minority = minority,
                                          MMProbabilityList = MMProbabilityList,
                                          higherClassColumn = higherClassColumn,
                                          crossValidation = crossValidation)

  return(predictionsMM)
}
