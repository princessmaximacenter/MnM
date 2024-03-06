extractIndividualValuesF1 <- function(predictionsMM,
                                      metaDataRef,
                                      classColumn,
                                      probabilityThreshold ,
                                      filterOrNot
                                      ) {


  nCases <- c(1,5,10,20,40,100)
  patientsPerTumor <- table(metaDataRef[,classColumn])

  if (filterOrNot == T) {
  predictionsMMFiltered <- predictionsMM %>% filter(probability1 > probabilityThreshold)
  } else {
    predictionsMMFiltered <- predictionsMM
  }
  tumorConfusionMatrix <- caret::confusionMatrix(factor(predictionsMMFiltered$predict,
                                                 levels = unique(c(predictionsMMFiltered$originalCall))),
                                          factor(predictionsMMFiltered$originalCall, levels =  unique(c(predictionsMMFiltered$originalCall))),
                                          dnn = c("Prediction", "Reference"))

  subsets <- c("3 <= n <= 5", paste( nCases[-c(length(nCases))],"< n <=",nCases[-c(1)])[-1], "n > 100")

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

    subsetResultTest <- predictionsMM %>% filter(originalCall %in% selectionMoreThanNCases)

    filteredResultTest <- subsetResultTest %>% filter(probability1 > probabilityThreshold)

    tumorConfusionMatrixSelection <- tumorConfusionMatrix$byClass %>% as.data.frame() %>% filter(rownames(.) %in% paste0("Class: ", selectionMoreThanNCases))

    tumorConfusionMatrixSelection$F1[is.na(tumorConfusionMatrixSelection$F1)] <- 0
    tumorConfusionMatrixSelection$Precision[is.na(tumorConfusionMatrixSelection$Precision)] <- 0
    recall <- table(filteredResultTest[,"originalCall"]) / patientsPerTumor[names(patientsPerTumor) %in% names(table(filteredResultTest[,"originalCall"]))]


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


  # # Now add the results if all samples are taken into account
  # #subsetResultTest <- predictionsMM %>% filter(originalCall %in% selectionMoreThanNCases)
  #
  #
  # tumorConfusionMatrixSelection <- tumorConfusionMatrix$byClass %>% as.data.frame() #%>% filter(rownames(.) %in% paste0("Class: ", selectionMoreThanNCases))
  #
  # tumorConfusionMatrixSelection$F1[is.na(tumorConfusionMatrixSelection$F1)] <- 0
  # tumorConfusionMatrixSelection$Precision[is.na(tumorConfusionMatrixSelection$Precision)] <- 0
  #
  # subsetResultTest <- predictionsMM
  #
  # filteredResultTest <- subsetResultTest %>% filter(probability1 > probabilityThreshold)
  #
  # recall <- table(filteredResultTest[,"originalCall"]) / patientsPerTumor[names(patientsPerTumor) %in% names(table(filteredResultTest[,"originalCall"]))]
  #
  # Precision <- tumorConfusionMatrixSelection[paste0("Class: ", names(recall)),"Precision"]
  # #Recall[i] <- mean(tumorConfusionMatrixSelection$Recall)
  # F1 <- tumorConfusionMatrixSelection[paste0("Class: ", names(recall)),"F1"]
  #
  # ourDF <- data.frame(nCases = "All",
  #                     tumorType = names(recall),
  #                     Precision = Precision,
  #                     F1 = F1,
  #                     Recall = as.numeric(recall))
  #
  #
  # totalDF <- rbind(totalDF, ourDF)



  # Summarize over 10 folds.


return(totalDF)
}
