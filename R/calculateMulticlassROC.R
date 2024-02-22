calculateMulticlassROC <- function(predictionsMMFinal,
                                   metaDataRef,
                                   classColumn,
                                   MMProbabilityList,
                                   filtering,
                                   probabilityScore = 0.8
                                   ) {

  #if (filtering == T) {
  #  predictionsMMFinal <- predictionsMMFinal %>% filter(probability1 > probabilityScore)
  #}
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
