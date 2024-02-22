shiftMedian <- function(countDataForScaling,
                        metaDataForScaling) {

  countDataForScaling %<>% as.data.frame()
  tumorTypes <- metaDataForScaling[,classColumn] %>% table() %>% names()
  for (i in seq(1:length(tumorTypes))) {
    metaDataSelected <- metaDataForScaling %>% filter(!!sym(classColumn) == tumorTypes[i], Status == "Primary")
    countDataSelected <- countDataForScaling[,rownames(metaDataSelected)]
    meanVals <- apply(countDataSelected, 1, mean)
    countDataSelected <- countDataSelected[meanVals >= 5,]
    countDataMedian <- apply(countDataSelected, 1, median)

    countDataScaled <- apply(countDataSelected, 2, function(x) x - countDataMedian) %>% as.data.frame()


    if (i == 1) {
      totalCountDataScaled <- countDataScaled
    } else {
      totalCountDataScaled <- totalCountDataScaled %>% filter(rownames(.) %in% rownames(countDataScaled))
      countDataScaled <- countDataScaled %>% filter(rownames(.) %in% rownames(totalCountDataScaled))
      totalCountDataScaled <- cbind(totalCountDataScaled,
                                    countDataScaled
                                    )
    }

  }




}
