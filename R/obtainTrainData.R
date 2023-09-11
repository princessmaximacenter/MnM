#' Obtain training data subsets
#'
#' This function is used to generate different training subsets for the majority voting system within the M&M setup.
#'
#' @param metaData Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param maxSamplesPerType How many samples should we maximally use per tumor (sub)type?
#' @param nModels How many models should be created for the majority voting system?
#'
#' @return A list with the specified number of different training data subsets from the reference cohort,
#' all with the specified maximum samples per tumor (sub)type.
#' @export
#'
obtainTrainData <- function(metaData, classColumn, maxSamplesPerType = 50, nModels = 100) {

  samplesTrainDefList <- list()
  metaData[, classColumn] <- as.factor(metaData[, classColumn])

  typesInFold <- table(metaData[, classColumn])
  typeProbs <- 1/sqrt((typesInFold))

  for ( j in c(1:nModels)){
    set.seed(j)
    samplesTrainDef <- c()

    samplesTrain <- unique(sample(rownames(metaData),
                                  prob=typeProbs[metaData[ , classColumn]],
                                  replace = T))

    includedTypes <- unique(metaData[samplesTrain, classColumn])

    # Look which tumor types are not present in your selected dataset
    missingType <- metaData[!(metaData[ , classColumn] %in% includedTypes), ]

    # If there are tumor types not present, we loop through them all to add
    # one entry for each tumor type to the training samples.
    if (length(rownames(missingType)) > 0) {
      missingTumors <- unique(missingType[,classColumn])

      for(i in seq(1:length(missingTumors))) {
        currentTumorType <- missingType[missingType[ , classColumn] == missingTumors[i],]
        samplesTrain <- c(samplesTrain, sample(rownames(currentTumorType), size = 1))
      }
    }

    #
    nSamplesPerType <- table(metaData[samplesTrain,classColumn])


    for (n in c(1:length(nSamplesPerType))){
      curSamples <- samplesTrain[metaData[samplesTrain,classColumn] ==
                                   names(nSamplesPerType)[n]]

      # Entries that are present more often than the specified maxSamplesPerType are downsampled to maxSamplesPerType
      if (nSamplesPerType[n] > maxSamplesPerType){

        samplesTrainDef <- c(samplesTrainDef,sample(curSamples,maxSamplesPerType,replace = F))

      }else{
        samplesTrainDef <- c(samplesTrainDef,curSamples)
      }
    }
    samplesTrainDefList[[j]] <- samplesTrainDef
  }
  return(samplesTrainDefList)
}
