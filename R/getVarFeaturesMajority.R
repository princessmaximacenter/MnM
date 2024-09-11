getVarFeaturesMajority <- function(
    dataTrain,
    samplesTrainDefList,
    nFeatures = 2500,
    nModels = 100

    ) {


  scaleFeaturesList <- base::list()
  for (i in base::seq(1:nModels)) {
    base::cat(base::paste("Working on model", i, "\n"))
    samplesTrainDef <- samplesTrainDefList[[i]]

    train.data <- dataTrain[, samplesTrainDef]

    varGenes <- base::apply(train.data,1 ,stats::var)

    varFeatures <- base::names(varGenes)[base::order(varGenes,decreasing = T)][c(1:nFeatures)]

    dataTrainFiltered <- base::as.data.frame(train.data) %>%
      dplyr::filter(base::rownames(.) %in% varFeatures)

    meanGenes <- base::apply(dataTrainFiltered, 1, base::mean,
                             na.rm = T)

    sdGenes <- base::apply(dataTrainFiltered, 1, stats::sd)

    scaleFeatures <- base::list(varFeatures = varFeatures,
                                meanGenes = meanGenes,
                                sdGenes = sdGenes)

    scaleFeaturesList[[i]] <- scaleFeatures

  }

  return(scaleFeaturesList)

}
