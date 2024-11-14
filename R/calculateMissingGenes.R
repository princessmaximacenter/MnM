#' Impute missing RNA-transcripts
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data for the reference cohort. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param countDataNew Matrix containing the RNA-transcript per million data for the new samples to be classified.
#' Samples are in the columns, different RNA-transcripts in the rows.
#' @param neededGenes RNA-transcripts that will be used within the classification process and are required for M&M to run.
#' Therefore, these genes will be imputed in case they are missing from countDataNew.
#' @param whichK The number of neighbor datapoints that need to be used to calculate missing RNA-transcripts from.
#'
#' @return Matrix containing the RNA-transcript per million data for the new samples to be classified,
#' now containing all the needed genes for the classification process.
calculateMissingGenes <- function(countDataRef,
                                  countDataNew,
                                  neededGenes,
                                  whichK = 3) {

  missingGenes <- neededGenes[neededGenes %notin% base::rownames(countDataNew)]

  missingDF <- matrix(NA, nrow = length(missingGenes), ncol = base::ncol(countDataNew))
  base::rownames(missingDF) <- missingGenes
  base::colnames(missingDF) <- base::colnames(countDataNew)
  newDataFrame <- base::rbind(countDataNew,
                        missingDF)

  countDataNewSubset <- newDataFrame[neededGenes,]
  countDataRefSubset <- countDataRef[neededGenes, ]

  hush=function(code){
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }

  # Subdivide dataset in chunks so that there are always more reference samples than samples to impute

  howManySplits <- base::ceiling(base::ncol(newDataFrame) / 500)

  splits <- base::rep(1:howManySplits, each = 500)
  splits <- splits[1:base::ncol(newDataFrame)]
  if (base::requireNamespace("impute") == F & base::requireNamespace("BiocManager") == F) {
    base::install.packages("BiocManager")
  } else if (base::requireNamespace("impute") == F) {
    BiocManager::install("impute")
  }

  for (i in unique(splits)) {
    combiSubset <- base::cbind(countDataRefSubset, countDataNewSubset[,splits == i])
    imputedDataComplete <- hush(impute::impute.knn(base::as.matrix(combiSubset), k = whichK))

    imputedData <- imputedDataComplete$data[missingGenes,base::colnames(countDataNewSubset[,splits == i])]

    if (i == 1) {
      imputedDataTotal <- imputedData
    } else {
      imputedDataTotal <- base::cbind(imputedDataTotal, imputedData)
    }
  }

  countDataNewTotal <- base::rbind(countDataNew,
                             imputedDataTotal)
  return(countDataNewTotal)
}
