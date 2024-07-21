#' Apply upsimpler on several sample sets, resampling them (for testing purposes)
#'
#' Function to apply the upsimpler algorithm on a list of sample sets.
#' Given a class of interest and a target sample from that class, this function
#' detects samples from this class and replaces them with new synthetic samples,
#' using the target sample as the seed. Whenever encountered, the target sample
#' is not replaced.
#'
#' @param sampleSets List of different training data sample subsets used for subsetting the available training dataset.
#' These contain the actual samples, not just the IDs (as they were selected by \code{selectFeaturesMajority}; for the
#' minority classifier, just use the reduced dataset in a single-element list).
#' @param metaDataRef Dataframe containing the metadata of the samples, with the sample IDs as rownames.
#' @param classColumn Column name in the metadata dataframe containing the class information.
#' @param sampleColumn Column name in the metadata dataframe containing the sample identifier.
#' @param targetSample A list with the target sample ID ($id) and actual data ($data) to be used as seed for the
#' upsimpler algorithm. The actual data are needed because some sample sets might not contain the target sample.
#' Make sure that the data are a dataframe, so that the features can be subsetted correctly.
#' @param upsimplerModule The upsimpler Python package as imported by \code{reticulate}.
#' @param upsimplerArgs List of arguments to be passed to the upsimpler algorithm, separately for the initialization of the
#' \code{Upsimpler} class ($init), as well as for the \code{upsimple} method ($upsimple).
#' @param whichSeed Seed to be used for the random number generator. Default is 1.
#'
#' @return A replicated list of the input sample sets with the samples from the class of interest replaced
#' by synthetic samples, except for the target sample. Note that the names of the synthetic samples are
#' replaced by the names produced by upsimpler.
#'
upsimplerResampling <- function(sampleSets,
                                metaDataRef,
                                classColumn,
                                sampleColumn,
                                targetSample,
                                upsimplerModule,
                                upsimplerArgs,
                                whichSeed = 1) {

  print <- function(x) base::print(paste('[APPLY_UPSIMPLER] ', x))

  targetSampleID <- targetSample$id
  targetClass <- metaDataRef[targetSampleID, classColumn]
  print(paste('Target sample:', targetSampleID, 'from class:', targetClass))

  for (i in base::seq_along(sampleSets)) {
    print(paste('Working on model', i))

    samples <- sampleSets[[i]]
    sampleIDs <- base::colnames(samples)
    # assert if there are sampleIDs that do NOT exist in the metadata
    if (any(!sampleIDs %in% base::rownames(metaDataRef))) {
      print(base::setdiff(sampleIDs, base::rownames(metaDataRef)))
      stop('The above sample IDs do not exist in the metadata.')
    }
    sampleClasses <- metaDataRef[sampleIDs, classColumn]
    nSamplesOfInterest <- sum(sampleClasses == targetClass)

    targetSampleIdx <- base::which(sampleIDs == targetSampleID)
    targetSampleFound <- length(targetSampleIdx) > 0
    if (targetSampleFound) {
      # we don't want to replace the target sample
      nSamplesOfInterest <- nSamplesOfInterest - 1
    } else {
      # we need to append the target sample to the samples
      thisSetFeatures <- base::rownames(samples)
      samples <- cbind(samples, targetSample$data[thisSetFeatures, , drop = F])
    }

    if (nSamplesOfInterest == 0) {
      print('  No samples of interest found in this sample set. Skipping upsimpler.')
      next
    }

    print(paste('  Found', nSamplesOfInterest, 'samples of interest in this sample set.'))
    if (targetSampleFound) print('  These include the target sample.')
    print(paste('  Replacing samples from class', targetClass, 'with synthetic samples...'))

    upsimplerArgs$init$dataDF <- samples %>% base::t() %>% base::as.data.frame()
    upsimplerArgs$init$metadataDF <- metaDataRef
    up <- do.call(upsimplerModule$Upsimpler, upsimplerArgs$init)

    upsimplerArgs$upsimple$n_samples <- as.integer(nSamplesOfInterest)
    upsimplerArgs$upsimple$rseed <- as.integer(whichSeed)
    synths <- do.call(up$upsimple, upsimplerArgs$upsimple) %>% base::t() %>% base::as.data.frame()

    # replace the samples from the class of interest with the synthetic samples
    # but do not replace the target sample
    # step 1: remove the samples from the class of interest but not the target sample
    sampleFilter <- sampleClasses == targetClass & sampleIDs != targetSampleID
    sampleSets[[i]] <- sampleSets[[i]][, !sampleFilter]
    # step 2: append the synthetic samples
    sampleSets[[i]] <- cbind(sampleSets[[i]], synths)
    # No need to deal with the metadata; from outside, it is very easily deduced
    print(paste('  Replacement of samples from class', targetClass, 'completed.'))
  }

  return(sampleSets)
}
