compareMethylationClassifier <- function(methylationClassifierData,
                                         predictionsMM,
                                         methylationThreshold = 0.9) {

  totalSamples <- nrow(methylationClassifierData)
  methClassifiedSamples <- rownames(methylationClassifierData)
  methErrors <- methylationClassifierData %>% filter(Disease_sub_class != methClass) %>% nrow()
  methCorrect <- methylationClassifierData %>% filter(Disease_sub_class == methClass) %>% nrow()
  methAccuracy <- round(methCorrect / totalSamples, digits = 2)

  methRecalled <- methylationClassifierData %>% filter(methClassScore > methylationThreshold) %>% nrow()
  methRecall <- methRecalled / totalSamples
  methPrecisioned <- methylationClassifierData %>% filter(methClassScore > methylationThreshold,
                                     Disease_sub_class == methClass) %>% nrow()
  methPrecision <- methPrecisioned / methRecalled

  MnMErrors <- predictionsMM %>% filter(rownames(.) %in% methClassifiedSamples,
                                        predict != originalCall) %>% nrow()

  MnMCorrect <- predictionsMM %>% filter(rownames(.) %in% methClassifiedSamples,
                                        predict == originalCall) %>% nrow()

  MnMAccuracy <- round(MnMCorrect / totalSamples, digits = 2)
  MnMRecalled <- predictionsMM %>% filter(probability1 > 0.8,
                                          rownames(.) %in% methClassifiedSamples) %>% nrow()
  MnMRecall <- MnMRecalled / totalSamples
  MnMPrecisioned <- predictionsMM %>% filter(probability1 > 0.8,
                                             predict == originalCall,
                                             rownames(.) %in% methClassifiedSamples) %>% nrow()
  MnMPrecision <- MnMPrecisioned / MnMRecalled

  methVSMnM <- data.frame(type = c("meth", "MnM"),
                          correct = c(methCorrect, MnMCorrect),
                          errors = c(methErrors, MnMErrors),
                          accuracy = c(methAccuracy, MnMAccuracy),
                          precision = c(methPrecision, MnMPrecision),
                          recall = c(methRecall, MnMRecall)
                          )
  return(methVSMnM)
}
