applyUpsimplerHeterogeneously <- function(samplesPerModel,
                                          metadataDF,
                                          classColumn,
                                          sampleColumn,
                                          seedSamples,
                                          upsimplerModule,
                                          upsimplerArgs) {
  synthsPerModel <- base::list()
  for (i in base::seq_along(samplesPerModel)) {
    upsimplerArgs$init$dataDF <- samplesPerModel[[i]] %>% base::t() %>% base::as.data.frame()
    upsimplerArgs$init$metadataDF <- metadataDF
    upsimplerArgs$init$class_col <- classColumn

    upsimpler <- do.call(upsimplerModule$Upsimpler, upsimplerArgs$init)

    synths <- applyUpsimpler(upsimpler = upsimpler,
                             upsimplerArgs = upsimplerArgs$upsimple,
                             targetSampleIDs = seedSamples,
                             metadataDF = metadataDF,
                             classColumn = classColumn,
                             sampleColumn = sampleColumn)

    # add the synths to the training data
    samplesPerModel[[i]] <- base::cbind(samplesPerModel[[i]], synths$synthDataDF %>% base::t() %>% base::as.data.frame())
    # append the synths metadata to the original metadata
    synths$synthMetadataDF[base::setdiff(base::names(metadataDF), base::names(synths$synthMetadataDF))] <- NA
    # append them row by row, only if the synths are not already present in the metadata
    for (j in base::seq_len(base::nrow(synths$synthMetadataDF))) {
      if (!(base::rownames(synths$synthMetadataDF)[j] %in% base::rownames(metadataDF))) {
        metadataDF <- base::rbind(metadataDF, synths$synthMetadataDF[j,])
      }
    }

    synthsPerModel[[i]] <- synths$synthDataDF
  }

  return(list(samplesPerModel = samplesPerModel,
              synthsPerModel = synthsPerModel,
              metadataDF = metadataDF))
}
