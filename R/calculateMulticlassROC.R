#' Calculate per tumor entity AUC curve
#'
#' @param predictionsMMFinal Dataframe containing the final classification ($predict) and original diagnosis label ($originalCall).
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param MMProbabilityList List with the probability scores for all classification labels that were assigned to samples.
#' @param filtering Do you want to use filtered classifications only?
#' @param probabilityThreshold What is the threshold you would like to use to call a classification 'confident'?
#'
#' @return Dataframe containing the per-sample probability scores for each tumor entity.
#'
calculateMulticlassROC <- function(predictionsMMFinal,
                                   metaDataRef,
                                   classColumn,
                                   MMProbabilityList,
                                   filtering,
                                   probabilityThreshold
                                   ) {

  if (filtering == T) {
   predictionsMMFinal <- predictionsMMFinal %>% filter(probability1 > probabilityThreshold)
  }
  probabilityForAUCFiltered <- data.frame(
    originalCall = predictionsMMFinal$originalCall)

  rownames(probabilityForAUCFiltered) <- rownames(predictionsMMFinal)

  #tumors <- unique(metaDataRef[,classColumn])
  tumors <- unique(predictionsMMFinal$originalCall)
  for (i in seq(1:length(tumors))) {
    probabilityForAUCFiltered[, tumors[i]] <- 0

  }

  if (filtering == T) {
  MMProbabilityListFiltered <- list()
  for (i in seq(1:length(MMProbabilityList))) {
    if (names(MMProbabilityList[i]) %in% rownames(predictionsMMFinal)) {
      MMProbabilityListFiltered[[names(MMProbabilityList[i])]] <- MMProbabilityList[[i]]
    }
  }
  } else {
    MMProbabilityListFiltered <- MMProbabilityList
  }
  for (i in seq(1:length(MMProbabilityListFiltered))) {
    predictions <- MMProbabilityListFiltered[[i]]
    patient <- names(MMProbabilityListFiltered[i])
    for (j in seq(1:length(predictions))) {
      if (names(predictions)[j] %in% tumors) {
      probabilityForAUCFiltered[patient, names(predictions)[j]] <- as.numeric(predictions[j])
      }
    }
  }

  if (filtering == T) {
  probabilityForAUCFiltered <- probabilityForAUCFiltered[, -c(56:60)]
  }

  "Make true diagram"
  trueDiagram <- matrix(0, nrow = nrow(probabilityForAUCFiltered),
                        ncol = ncol(probabilityForAUCFiltered) - 1) %>% as.data.frame()

  colnames(trueDiagram) <- colnames(probabilityForAUCFiltered)[-c(1)]
  rownames(trueDiagram) <- rownames(probabilityForAUCFiltered)

  for (i in seq(1:nrow(trueDiagram))) {
    originalCall <- probabilityForAUCFiltered$originalCall[i]
    trueDiagram[i,originalCall] <- 1
  }

  colnames(trueDiagram) <- paste0(colnames(probabilityForAUCFiltered)[-c(1)], "_true")
  probabilityForAUC_PR <- probabilityForAUCFiltered
  colnames(probabilityForAUC_PR)[-c(1)] <- paste0(colnames(probabilityForAUCFiltered)[-c(1)], "_pred_RF")

  colnames(probabilityForAUC_PR)[-c(1)] <- gsub(" ", ".", colnames(probabilityForAUC_PR)[-c(1)])
  colnames(trueDiagram) <- gsub(" ", ".", colnames(trueDiagram))

  prDFFiltered <- cbind(trueDiagram, probabilityForAUC_PR)

  return(prDFFiltered)


}
