#' Calculate per tumor (sub)type performance
#'
#' @param predictionsMM Predictions for samples by M&M, coming from the function integrateMM.
#' @param metaDataRef  Metadata file containing the links between the samples and
#' the tumor (sub)type diagnosis within the reference cohort.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param probabilityThreshold What is the probability score threshold you would like to use to call a classification 'confident'?
#' @param filterOrNot Do you want to filter the 'confident' classifications only for your calculation?
#' @import caret
#'
#' @return Dataframe containing the tumor (sub)type ($tumorType) with its associated precision ($Precision), F1 score ($F1),
#' sensitivity ($Sensitivity) and recall ($Recall), stratified by population frequency ($nCases).
#'
#' Note that recall here means the fraction of samples obtaining a 'confident' classification
#' (probability higher than probabilityThreshold).
#'

extractIndividualValuesF1 <- function(predictionsMM,
                                      metaDataRef,
                                      classColumn,
                                      probabilityThreshold,
                                      filterOrNot
                                      ) {


  nCases <- c(1,5,10,20,40,100)
  patientsPerTumor <- base::table(metaDataRef[,classColumn])

  if (filterOrNot == T) {
    predictionsMMFiltered <- predictionsMM %>% dplyr::filter(probability1 > probabilityThreshold)
  } else {
    predictionsMMFiltered <- predictionsMM
  }

  tumorConfusionMatrix <- caret::confusionMatrix(factor(predictionsMMFiltered$predict,
                                                 levels = base::unique(c(predictionsMMFiltered$originalCall))),
                                          factor(predictionsMMFiltered$originalCall,
                                                 levels =  base::unique(c(predictionsMMFiltered$originalCall))),
                                          dnn = c("Prediction", "Reference"))

  subsets <- c("3 <= n <= 5", paste( nCases[-c(base::length(nCases))],"< n <=",nCases[-c(1)])[-1], "n > 100")

  for ( i in 1:(length(nCases))){

    if ( i == length(nCases)){
      # Select tumors within a certain block (1-3, 4-5, 6-10, etc.)
      moreThanNCases <-  patientsPerTumor > nCases[i]
    } else {
      moreThanNCases <- (patientsPerTumor > nCases[i] &
                           patientsPerTumor <= nCases[i+1])
    }
    # Find their names

    selectionMoreThanNCases <- names(patientsPerTumor[moreThanNCases])

    subsetResultTest <- predictionsMM %>% dplyr::filter(originalCall %in% selectionMoreThanNCases)

    filteredResultTest <- subsetResultTest %>% dplyr::filter(probability1 > probabilityThreshold)

    tumorConfusionMatrixSelection <- tumorConfusionMatrix$byClass %>%
      as.data.frame() %>%
      dplyr::filter(rownames(.) %in% paste0("Class: ", selectionMoreThanNCases))

    tumorConfusionMatrixSelection$F1[is.na(tumorConfusionMatrixSelection$F1)] <- 0
    tumorConfusionMatrixSelection$Precision[is.na(tumorConfusionMatrixSelection$Precision)] <- 0
    filteredResultTest[,"originalCall"] <- factor(filteredResultTest[,"originalCall"], levels = selectionMoreThanNCases)
    recall <- base::table(filteredResultTest[,"originalCall"]) /
      patientsPerTumor[base::names(patientsPerTumor) %in% base::names(base::table(filteredResultTest[,"originalCall"]))]


    Precision <- tumorConfusionMatrixSelection[paste0("Class: ", names(recall)),"Precision"]
    #Recall[i] <- mean(tumorConfusionMatrixSelection$Recall)
    F1 <- tumorConfusionMatrixSelection[paste0("Class: ", names(recall)),"F1"]
    Sensitivity <- tumorConfusionMatrixSelection[paste0("Class: ", names(recall)),"Sensitivity"]
    ourDF <- data.frame(nCases = subsets[i],
                        tumorType = names(recall),
                        Precision = Precision,
                        F1 = F1,
                        Recall = as.numeric(recall),
                        Sensitivity = Sensitivity
                        )

    if (i == 1) {
      totalDF <- ourDF
    } else {
      totalDF <- rbind(totalDF, ourDF)
    }
  }

return(totalDF)
}
