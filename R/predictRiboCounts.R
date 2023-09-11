#' Predict the
#'
#' @param riboModel
#' @param data
#'
#' @return
#' @export
#'
predictRiboCounts <- function(riboModel, data) {
  # Predict how much protein coding reads v.s. the ribosomal reads are present within the data
  relevantCoefficients <- riboModelCounts$relevantCoefficients
  highMeanGenes <- names(riboModelCounts$meanGenes)
  meanGenes <- riboModelCounts$meanGenes
  varGenes <- riboModelCounts$varGenes
  normalizedData <- apply(data,2,function(x) (x[highMeanGenes]-meanGenes[highMeanGenes])/varGenes[highMeanGenes])

  predictProteinCoding <- 1-(relevantCoefficients[1] + t(normalizedData[names(relevantCoefficients)[-1],]) %*% relevantCoefficients[-1])

  # Divide counts for each sample by the percentage of protein coding sequence.
  # The less protein coding reads, the higher the scale-up of eventual reads.
  data <- sapply(c(1:ncol(data)),function(x) smart.round(data[,x]/predictProteinCoding[x],digits =3))
  colnames(data) <- colnames(normalizedData)

  # After scaling-up samples with limited numbers of protein coding reads,
  # the most important ribosomal proteins are removed from the dataset.
  data <- data[!(rownames(data) %in% names(relevantCoefficients)[-1]),]
  return(data)
}
