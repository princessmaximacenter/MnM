#' Obtain training data subsets
#'
#' This function is used to generate different training subsets for the different models within the M&M setup.
#'
#' @param metaDataRef Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param maxSamplesPerType How many samples should we maximally use per tumor subtype?
#' @param nModels How many models should be created for the classifier?
#'
#' @return A list with the specified number of different training data subsets from the available complete training data,
#' all with the specified maximum samples per tumor subtype.
#'
obtainTrainData <- function(metaDataRef, classColumn, maxSamplesPerType = 50, nModels = 100) {

  samplesTrainDefList <- list()
  metaDataRef[, classColumn] <- as.factor(metaDataRef[, classColumn])

  typesInFold <- base::table(metaDataRef[, classColumn])
  typeProbs <- 1/sqrt(typesInFold)

  for ( j in c(1:nModels)){
    #set.seed(j)
    samplesTrainDef <- c()

    samplesTrain <- base::unique(base::sample(rownames(metaDataRef),
                                  prob=typeProbs[metaDataRef[ , classColumn]],
                                  replace = T))

    includedTypes <- base::unique(metaDataRef[samplesTrain, classColumn])

    # Look which tumor types are not present in your selected dataset
    missingType <- metaDataRef[!(metaDataRef[ , classColumn] %in% includedTypes), ]

    # If there are tumor types not present, we loop through them all to add
    # one entry for each tumor type to the training samples.
    if (base::length(rownames(missingType)) > 0) {
      missingTumors <- base::unique(missingType[,classColumn])

      for(i in seq(1:base::length(missingTumors))) {
        currentTumorType <- missingType[missingType[ , classColumn] == missingTumors[i],]
        samplesTrain <- c(samplesTrain, base::sample(rownames(currentTumorType), size = 1))
      }
    }

    #
    nSamplesPerType <- base::table(metaDataRef[samplesTrain,classColumn])


    for (n in c(1:base::length(nSamplesPerType))){
      curSamples <- samplesTrain[metaDataRef[samplesTrain,classColumn] ==
                                   base::names(nSamplesPerType)[n]]

      # Entries that are present more often than the specified maxSamplesPerType are downsampled to maxSamplesPerType
      if (nSamplesPerType[n] > maxSamplesPerType){

        samplesTrainDef <- c(samplesTrainDef,base::sample(curSamples,maxSamplesPerType,replace = F))

      }else{
        samplesTrainDef <- c(samplesTrainDef,curSamples)
      }
    }
    samplesTrainDefList[[j]] <- samplesTrainDef
  }
  return(samplesTrainDefList)
}
