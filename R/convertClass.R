convertClass <- function(classifierResults,
                         classColumn,
                         higherClassColumn,
                         metaDataRef,
                         nModels = 100,
                         crossValidation = F) {


 probabilities <- obtainProbabilities(classifierResults,
                                         nModels = nModels,
                                         crossValidation = crossValidation)
  probabilityDF <- obtainProbabilityDF(classifierResults,
                                       probabilities)

  probabilityDFWithHigherClass <- findPredictedHigherClass(probabilityDF,
                                                           metaDataRef,
                                                           classColumn,
                                                           higherClassColumn,
                                                           crossValidation = crossValidation)

  higherClassList <- obtainOtherScores(classifierResults,
                                       metaDataRef,
                                       classColumn = classColumn,
                                       higherClassColumn = higherClassColumn,
                                       nModels = nModels,
                                       crossValidation = crossValidation)

  resultHigherClass <- addHigherClassScore(probabilityDFWithHigherClass, higherClassList)
  resultHigherClass <- integrateScores(resultHigherClass)
  return(resultHigherClass)
}



checkPredictions <- function(classifications, metaDataTest, classColumn, higherClassColumn, metaDataRef) {
  # Look at the original calls for each test sample
  originalCall <- metaDataTest[rownames(classifications),classColumn]
  originalCallHigherClass <- metaDataTest[rownames(classifications),higherClassColumn]
  # Store the bestFit, the Originalcall and the accompanying probability score within the final dataframe.
  predictionWithOriginalCalls <- cbind(classifications,
                                       originalCall = originalCall,
                                       originalCallHigherClass = originalCallHigherClass)




  return(predictionWithOriginalCalls)
  # Make sure that the ultimatePredictions still have their accompanying biomaterial_i
}
