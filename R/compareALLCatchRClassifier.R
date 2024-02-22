compareALLCatchRClassifier <- function(predictionsMM,
                                       predictionsAllCatchR,
                                       metaDataRef,
                                       classColumn,
                                       MnMThreshold = 0.7,
                                       whichSeed
                                      # howManyFolds = 1
                                      ) {


  set.seed(whichSeed)
  predictionsAllCatchRTotal <- predictionsAllCatchR
  #for ( i in seq(1:howManyFolds)) {
  BALL <- metaDataRef[rownames(predictionsAllCatchRTotal),] %>% filter(Disease_sub_class == "B-ALL",
                                 Disease_sub_specification1 != "NOS B-ALL") %>% rownames(.)
  #if (howManyFolds != 1) {
  #BALL <- unique(sample(BALL, replace = T))
  #}
  #BALL <- sample(BALL, size = 100)
  totalSamples <- length(BALL)
  predictionsAllCatchR <- predictionsAllCatchRTotal %>% filter(rownames(.) %in% BALL)


  predictionsAllCatchR$originalCall <- metaDataRef[rownames(predictionsAllCatchR),classColumn ]


  AllCatchRCorrect <- predictionsAllCatchR %>% filter(originalCall == Prediction) %>% nrow()
  AllCatchRErrors <- predictionsAllCatchR %>% filter(originalCall != Prediction) %>% nrow()
  AllCatchRAccuracy <-  round(AllCatchRCorrect / totalSamples, digits = 2)
  AllCatchRRecalled <- predictionsAllCatchR %>% filter(Confidence == "high-confidence") %>% nrow()
  AllCatchRRecall <- AllCatchRRecalled / totalSamples
  AllCatchRPrecisioned <- predictionsAllCatchR %>% filter(Confidence == "high-confidence",
                                                          originalCall == Prediction) %>% nrow()
  AllCatchRPrecision <- AllCatchRPrecisioned / AllCatchRRecalled


  MnMCorrect <- predictionsMM %>% filter(rownames(.) %in% BALL,
                                         originalCall == predict) %>% nrow()
  MnMErrors <- predictionsMM %>% filter(rownames(.) %in% BALL,
                                  originalCall != predict) %>% nrow()
  MnMAccuracy <- round(MnMCorrect / totalSamples, digits = 2)
  MnMRecalled <- predictionsMM %>% filter(probability1 > MnMThreshold,
                                          rownames(.) %in% BALL) %>% nrow()

  MnMRecall <- MnMRecalled / totalSamples
  MnMPrecisioned <- predictionsMM %>% filter(probability1 > MnMThreshold,
                                             predict == originalCall,
                                             rownames(.) %in% BALL) %>% nrow()
  MnMPrecision <- MnMPrecisioned / MnMRecalled


  AllCatchRVSMnM <- data.frame(type = c("AllCatchR", "MnM"),
                          correct = c(AllCatchRCorrect, MnMCorrect),
                          errors = c(AllCatchRErrors, MnMErrors),
                          accuracy = c(AllCatchRAccuracy, MnMAccuracy),
                          precision = c(AllCatchRPrecision, MnMPrecision),
                          recall = c(AllCatchRRecall, MnMRecall)
  )
 # if (i == 1) {
 #   AllCatchRVSMnMTotal <- AllCatchRVSMnM
 # } else {
 #   AllCatchRVSMnMTotal <- rbind(AllCatchRVSMnMTotal, AllCatchRVSMnM)
 # }
#}
  return(AllCatchRVSMnM)
}
