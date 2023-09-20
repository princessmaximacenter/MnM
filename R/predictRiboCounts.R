#' Correct new count data for ribodepletion efficacy
#'
#' This function is designed to use the same ribodepletion correcting model
#' on the new samples as was used on the reference cohort.
#'
#' @param riboModel Model for ribocorrection, normally stored inside the R-object
#' obtained form running the function createScalingsMajority.
#' @param data Count data to be transformed.
#'
#' @return Count data corrected for the efficacy of the ribodepletion protocol.
#'
#'
predictRiboCounts <- function(riboModel, data) {
  # Predict how much protein coding reads v.s. the ribosomal reads are present within the data
  relevantCoefficients <- riboModel$relevantCoefficients
  highMeanGenes <- names(riboModel$meanGenes)
  meanGenes <- riboModel$meanGenes
  varGenes <- riboModel$varGenes
  dataSub <- data[highMeanGenes, , drop = F]
  normalizedData <- apply(dataSub,2,function(x) (x -meanGenes[highMeanGenes])/varGenes[highMeanGenes])

  predictProteinCoding <- 1-(relevantCoefficients[1] + t(normalizedData[names(relevantCoefficients)[-1], , drop = F]) %*% relevantCoefficients[-1])

  # Divide counts for each sample by the percentage of protein coding sequence.
  # The less protein coding reads, the higher the scale-up of eventual reads.

  data <- sapply(c(1:ncol(data)),function(x) smartRound(data[,x]/predictProteinCoding[x],digits =3))
  colnames(data) <- colnames(normalizedData)

  # After scaling-up samples with limited numbers of protein coding reads,
  # the most important ribosomal proteins are removed from the dataset.
  data <- data[!(rownames(data) %in% names(relevantCoefficients)[-1]), , drop = F]
  return(data)
}
