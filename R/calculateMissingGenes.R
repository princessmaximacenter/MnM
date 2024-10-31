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

  hush=function(code){
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }

  # Subdivide dataset in chunks so that there are always more reference samples than samples to impute

  howManySplits <- ceiling(ncol(newDataFrame) / 500)

  splits <- rep(1:howManySplits, each = 500)
  splits <- splits[1:ncol(newDataFrame)]

  for (i in unique(splits)) {
    combiSubset <- cbind(countDataRefSubset, countDataNewSubset[,splits == i])
    imputedDataComplete <- hush(impute::impute.knn(as.matrix(combiSubset), k = whichK))

    imputedData <- imputedDataComplete$data[missingGenes,colnames(countDataNewSubset[,splits == i])]

    if (i == 1) {
      imputedDataTotal <- imputedData
    } else {
      imputedDataTotal <- cbind(imputedDataTotal, imputedData)
    }
  }

  countDataNewTotal <- rbind(countDataNew,
                             imputedDataTotal)
  return(countDataNewTotal)
}
