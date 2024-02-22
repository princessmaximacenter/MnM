newBarplot <- function(metaDataRef,
                       classColumn,
                       predictionsMM,
                       probabilityThreshold) {


  nCases <- c(1,3,5,10,20,40,100)
  patientsPerTumor <- table(metaDataRef[,classColumn])
  nCasesText <- c("n = 3", paste( nCases[-length(nCases)],"< n <=",nCases[-1])[-1], "n > 100")

  # Need to add the tumors that are not in the metadata, but are in the test set

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

    for (j in seq(1:length(selectionMoreThanNCases))) {
      tumorType <- selectionMoreThanNCases[j]
      subsetResultTest <- predictionsMM %>% filter(originalCall == tumorType)
      filteredResultTest <- subsetResultTest %>% filter(probability1 > probabilityThreshold)
      tumorDF <- data.frame(tumorType = tumorType,
                            fractionClassified = round(nrow(filteredResultTest) / nrow(subsetResultTest), digits = 2),
                            nCases = nCasesText[i])

      if (i == 1 & j == 1) {
        tumorDFTotal <- tumorDF
      } else {
        tumorDFTotal <- rbind(tumorDFTotal,
                              tumorDF)
      }
    }

  }
  tumorDFTotal$nCases <- factor(tumorDFTotal$nCases, levels = unique(tumorDFTotal$nCases))
  tumorDFTotal %<>% arrange(nCases, fractionClassified)
  tumorDFTotal$tumorType <- factor(tumorDFTotal$tumorType,
                                         levels = unique(tumorDFTotal$tumorType))
  return(tumorDFTotal)
}
