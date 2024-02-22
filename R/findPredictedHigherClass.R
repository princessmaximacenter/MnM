findPredictedHigherClass <- function(probabilityDF, metaData, classColumn, higherClassColumn, crossValidation = T) {
  classifications <- probabilityDF

  linkClassAndHigherClass <- metaData[ , c(classColumn, higherClassColumn)] %>% unique

  predictWithHigherClass <- left_join(classifications, linkClassAndHigherClass, by = c("predict" = classColumn))


  colnumber <- which(colnames(predictWithHigherClass) == higherClassColumn)
  colnames(predictWithHigherClass)[colnumber] <- "predictHigherClass"

  predictWithHigherClass <- left_join(predictWithHigherClass, linkClassAndHigherClass, by = c("predict2" = classColumn))
  colnumber <- which(colnames(predictWithHigherClass) == higherClassColumn)
  colnames(predictWithHigherClass)[colnumber] <- "predictHigherClass2"

  predictWithHigherClass <- left_join(predictWithHigherClass, linkClassAndHigherClass, by = c("predict3" = classColumn))
  colnumber <- which(colnames(predictWithHigherClass) == higherClassColumn)
  colnames(predictWithHigherClass)[colnumber] <- "predictHigherClass3"

  rownames(predictWithHigherClass) <- rownames(classifications)
  if (crossValidation == T) {
    probabilityDFWithHigherClass <- left_join(predictWithHigherClass, linkClassAndHigherClass, by = c("originalCall" = classColumn))

    colnumber <- which(colnames(probabilityDFWithHigherClass) == higherClassColumn)
    colnames(probabilityDFWithHigherClass)[colnumber] <- "originalCallHigherClass"

    rownames(probabilityDFWithHigherClass) <- rownames(predictWithHigherClass)
  } else {
    probabilityDFWithHigherClass <- predictWithHigherClass
  }

  if (crossValidation == T) {
    probabilityDFWithHigherClass <- probabilityDFWithHigherClass[,c("originalCall",
                                                                    "predict",
                                                                    "probability1",
                                                                    "predict2",
                                                                    "probability2",
                                                                    "predict3",
                                                                    "probability3",
                                                                    "originalCallHigherClass",
                                                                    "predictHigherClass",
                                                                    "predictHigherClass2",
                                                                    "predictHigherClass3")]
  } else {
    probabilityDFWithHigherClass <- probabilityDFWithHigherClass[,c("predict",
                                                                    "probability1",
                                                                    "predict2",
                                                                    "probability2",
                                                                    "predict3",
                                                                    "probability3",
                                                                    "predictHigherClass",
                                                                    "predictHigherClass2",
                                                                    "predictHigherClass3")]
  }
  return(probabilityDFWithHigherClass)
}
