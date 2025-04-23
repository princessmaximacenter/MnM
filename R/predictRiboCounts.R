#' Correct new count data for ribodepletion efficacy
#'
#' This function is designed to use the same ribodepletion correcting model
#' on the new samples as was used on the reference cohort.
#'
#' @param riboModel Model for ribocorrection, normally stored inside the R-object
#' obtained from running the function createModelsMinority or createScalingsMajority as 'riboModelList'.
#' @param data Count data to be transformed.
#' @param countDataRef Matrix containing the RNA-transcript per million data. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param whichKimputation The number of neighbor datapoints that need to be used to calculate missing RNA-transcripts from (in case needed).
#' @return Count data corrected for the efficacy of the ribodepletion protocol.
#' @export
#'
predictRiboCounts <- function(riboModel, data, countDataRef, whichKimputation) {
  # Predict how much protein coding reads v.s. the ribosomal reads are present within the data
  `%notin%`  <- base::Negate(`%in%`)
  relevantCoefficients <- riboModel$relevantCoefficients
  highMeanGenes <- base::names(riboModel$meanGenes[names(relevantCoefficients[-1])])
  meanGenes <- riboModel$meanGenes[names(relevantCoefficients[-1])]
  varGenes <- riboModel$varGenes[names(relevantCoefficients[-1])]
  dataSub <- data[highMeanGenes[highMeanGenes %in% rownames(data)], , drop = F]
  if (sum(names(varGenes) %notin% rownames(data)) > 0) {

    missingGenes <- names(varGenes)[names(varGenes) %notin% rownames(dataSub)]
    cat(paste0("There are ", length(missingGenes), " genes missing from the dataset for ribodepletion correction.\nImputing their values.\n"))

    dataSub <- calculateMissingGenes(countDataNew = dataSub,
                                     neededGenes = names(varGenes),
                                     countDataRef = countDataRef,
                                     whichK = whichKimputation)
  }

  dataSub <- dataSub[highMeanGenes, , drop = F]
  normalizedData <- base::apply(dataSub,2,function(x) (x -meanGenes[highMeanGenes])/varGenes[highMeanGenes])

  predictProteinCoding <- 1-(relevantCoefficients[1] + base::t(normalizedData) %*% relevantCoefficients[-1])

  if (length(which(predictProteinCoding < 0)) > 0) {
    data <- data[,-c(which(predictProteinCoding < 0))]
    normalizedData <- normalizedData[,-c(which(predictProteinCoding < 0))]
    print(paste0("Found ", length(which(predictProteinCoding < 0)), " bad quality samples where there appears to be no protein content. Removing said samples from the classification procedure."))
    predictProteinCoding <- predictProteinCoding[colnames(data), , drop = F]
  }
  # Divide counts for each sample by the percentage of protein coding sequence.
  # The less protein coding reads, the higher the scale-up of eventual reads.

  data <- base::sapply(c(1:base::ncol(data)),function(x) smartRound(data[,x]/predictProteinCoding[x],digits =3))
  base::colnames(data) <- base::colnames(normalizedData)

  # After scaling-up samples with limited numbers of protein coding reads,
  # the most important ribosomal RNAs are removed from the dataset.
  data <- data[!(base::rownames(data) %in% base::names(relevantCoefficients)[-1]), , drop = F]
  return(data)
}
