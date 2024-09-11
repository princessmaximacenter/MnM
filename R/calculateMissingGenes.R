#' Impute missing RNA-transcripts
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data for the reference cohort. Samples are in the columns,
#' different RNA-transcripts in the rows.
#' @param countDataNew Matrix containing the RNA-transcript per million data for the new samples to be classified.
#' Samples are in the columns, different RNA-transcripts in the rows.
#' @param neededGenes RNA-transcripts that will be used within the classification process and are required for M&M to run.
#' Therefore, these genes will be imputed in case they are missing from countDataNew.
#'
#' @return Matrix containing the RNA-transcript per million data for the new samples to be classified,
#' now containing all the needed genes for the classification process.
calculateMissingGenes <- function(countDataRef,
                                  countDataNew,
                                  neededGenes,
                                  whichK = 10) {

  missingGenes <- neededGenes[neededGenes %notin% rownames(countDataNew)]

  missingDF <- matrix(NA, nrow = length(missingGenes), ncol = ncol(countDataNew))
  rownames(missingDF) <- missingGenes
  colnames(missingDF) <- colnames(countDataNew)
  newDataFrame <- rbind(countDataNew,
                        missingDF)

  countDataNewSubset <- newDataFrame[neededGenes,]
  countDataRefSubset <- countDataRef[neededGenes, ]

  combiSubset <- cbind(countDataRefSubset, countDataNewSubset)

  hush=function(code){
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }
  imputedDataComplete <- hush(impute::impute.knn(as.matrix(combiSubset), k = whichK))

  imputedData <- imputedDataComplete$data[missingGenes,colnames(countDataNewSubset)]
  countDataNewTotal <- rbind(countDataNew,
                              imputedData)
  return(countDataNewTotal)
}
