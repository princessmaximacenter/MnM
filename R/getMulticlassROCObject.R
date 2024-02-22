getMulticlassROCObject <- function(metaDataRef,
         classColumn,
         MMProbabilityList
) {

  probabilityPerClass <- data.frame(patient = rownames(metaDataRef),
                      originalCall = metaDataRef[, classColumn])

  rownames(probabilityPerClass) <- probabilityPerClass$patient
  tumors <- unique(metaDataRef$Disease_sub_class)

  for (i in seq(1:length(tumors))) {
    probabilityPerClass[, tumors[i]] <- 0

  }

  for (i in seq(1:length(MMProbabilityList))) {
    predictions <- MMProbabilityList[[i]]
    patient <- names(MMProbabilityList[i])
    for (j in seq(1:length(predictions))) {
      probabilityPerClass[patient, names(predictions)[j]] <- as.numeric(predictions[j])
    }
  }

  return(probabilityPerClass)

}


