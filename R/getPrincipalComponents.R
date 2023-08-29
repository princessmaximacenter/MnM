getPrincipalComponents <- function(dataTrain,
         metaData,
         samplesTrainDefList,
         classColumn,
         nModels = 100,
         nFeatures = 2500,
         nComps = 100
) {

  # The genes are in the columns right now, so the variance should be determined per column (removing)

  prList <- list()
  scaleFeaturesList <- list()
  for (i in seq(1:nModels)) {
    print(paste("Working on model", i))
    samplesTrainDef <- samplesTrainDefList[[i]]

    train.data <- dataTrain[, samplesTrainDef]

    varGenes <- apply(train.data,1 ,var)

    varFeatures <- names(varGenes)[order(varGenes,decreasing = T)][c(1:nFeatures)]

    dataTrainFiltered <- as.data.frame(train.data) %>%
      filter(rownames(.) %in% varFeatures)

    meanGenes <- apply(dataTrainFiltered, 1, mean,
                       na.rm = T)


    sdGenes <- apply(dataTrainFiltered, 1, sd)

    dataScale <- apply(dataTrainFiltered[varFeatures,], 2, function(x) (x-meanGenes[varFeatures])/sdGenes[varFeatures])

    pr <- prcomp(t(dataScale),
                 rank.=nComps)

    prList[[i]] <- pr
    scaleFeatures <- list(varFeatures = varFeatures,
                          meanGenes = meanGenes,
                          sdGenes = sdGenes)
    scaleFeaturesList[[i]] <- scaleFeatures

    rotationsAndScalingsList <- list(prList = prList,
                                     scaleFeaturesList = scaleFeaturesList)
  }
  return(rotationsAndScalingsList)
}
