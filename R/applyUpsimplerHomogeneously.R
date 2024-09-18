applyUpsimplerHomogeneously <- function(dataDF,
                                        metadataDF,
                                        samplesPerModel,
                                        classColumn,
                                        sampleColumn,
                                        seedSamples,
                                        nModels,
                                        maxSamplesPerType,
                                        upsimplerModule,
                                        upsimplerArgs) {

  allMajorityFeatures <- base::lapply(samplesPerModel, base::rownames) %>% base::unlist() %>% base::unique()

  # Apply the upsimpler algorithm
  upsimplerArgs$init$dataDF <- dataDF[allMajorityFeatures,] %>% base::t() %>% base::as.data.frame()
  upsimplerArgs$init$metadataDF <- metadataDF
  upsimplerArgs$init$class_col <- classColumn

  upsimpler <- do.call(upsimplerModule$Upsimpler, upsimplerArgs$init)

  synths <- applyUpsimpler(upsimpler = upsimpler,
                           upsimplerArgs = upsimplerArgs$upsimple,
                           targetSampleIDs = seedSamples,
                           metadataDF = metadataDF,
                           classColumn = classColumn,
                           sampleColumn = sampleColumn)

  # distribute the synths to the models
  synthsPerModel <- obtainTrainData(metaDataRef = synths$synthMetadataDF,
                                    classColumn = classColumn,
                                    nModels = nModels,
                                    maxSamplesPerType = maxSamplesPerType)

  # map synthsPerModel from sample IDs to the actual data
  synths$synthDataDF %<>% base::t() %>% base::as.data.frame()
  for (i in base::seq_along(synthsPerModel)) {
    synthsPerModel[[i]] <- synths$synthDataDF[base::rownames(samplesPerModel[[i]]), synthsPerModel[[i]]]
  }
  synths$synthDataDF %<>% base::t() %>% base::as.data.frame()

  # add the synths to the training data
  samplesPerModel <- base::lapply(base::seq_along(samplesPerModel), function(i) {
    base::cbind(samplesPerModel[[i]], synthsPerModel[[i]])
  })

  # synths metadata is a single column; we need to add the other columns
  # in order to rbind it to the original metadata
  synths$synthMetadataDF[base::setdiff(base::names(metadataDF), base::names(synths$synthMetadataDF))] <- NA
  metadataDF <- base::rbind(metadataDF, synths$synthMetadataDF)

  return(list(samplesPerModel = samplesPerModel,
              synthsPerModel = synthsPerModel,
              metadataDF = metadataDF))
}
