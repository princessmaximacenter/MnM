#' Correct new count data for ribodepletion efficacy
#'
#' This function is designed to use the same ribodepletion correcting model
#' on the new samples as was used on the reference cohort.
#'
#' @param riboModel Model for ribocorrection, normally stored inside the R-object
#' obtained from running the function createModelsMinority or createScalingsMajority as 'riboModelList'.
#' @param data Count data to be transformed.
#'
#' @return Count data corrected for the efficacy of the ribodepletion protocol.
#'
#'
predictRiboCounts <- function(riboModel, data) {
  # Predict how much protein coding reads v.s. the ribosomal reads are present within the data
  relevantCoefficients <- riboModel$relevantCoefficients
  highMeanGenes <- base::names(riboModel$meanGenes)
  meanGenes <- riboModel$meanGenes
  varGenes <- riboModel$varGenes
  dataSub <- data[highMeanGenes, , drop = F]
  normalizedData <- base::apply(dataSub,2,function(x) (x -meanGenes[highMeanGenes])/varGenes[highMeanGenes])

  predictProteinCoding <- 1-(relevantCoefficients[1] + base::t(normalizedData[base::names(relevantCoefficients)[-1], , drop = F]) %*% relevantCoefficients[-1])

  # Divide counts for each sample by the percentage of protein coding sequence.
  # The less protein coding reads, the higher the scale-up of eventual reads.

  data <- base::sapply(c(1:base::ncol(data)),function(x) smartRound(data[,x]/predictProteinCoding[x],digits =3))
  base::colnames(data) <- base::colnames(normalizedData)

  # After scaling-up samples with limited numbers of protein coding reads,
  # the most important ribosomal RNAs are removed from the dataset.
  data <- data[!(base::rownames(data) %in% base::names(relevantCoefficients)[-1]), , drop = F]
  return(data)
}
