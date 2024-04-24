#' Combine individual model results into final classification
#'
#' Convert individual classifications to an object containing the final per-sample classifications,
#' per-label probabilities and metadata that the classifications were based on.
#'
#' @param result Dataframe containing the classification labels from each model for each sample.
#' @param metaDataRef Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param addOriginalCall Would you like to add the original tumor subtype label from the metadata (TRUE)?
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#'
#' @return List containing the final per-sample classifications ($classifications),
#' per-label probabilities ($probabilityList) and metadata that the classifications were based on ($metaDataRef).

convertResultToClassification <- function(result,
         metaDataRef,
         addOriginalCall,
         classColumn = NA) {

  if (nrow(result) == 1 || typeof(apply(result, 1, table)) == 'integer') {
    randomVector <- paste0("fake", 1:ncol(result)) %>% as.data.frame() %>% t() %>% as.data.frame()
    colnames(randomVector) <- colnames(result)
    result1 <- rbind(result, randomVector)
    probability <-  apply(result1, 1, table)
    probability <- probability[-base::length(probability)]
  } else {
    # Find out how often a certain tumor type prediction is made for a specific sample
    probability <- apply(result, 1, base::table)
  }
  # Locate the position of the highest probability
  positions <- lapply(probability, which.max)
  positions <- unlist(positions)

  bestFit <- data.frame(predict = rep(NA, times = base::length(result$fold1)))
  probabilityScores <- vector()

  # Extract the different calls being made for each sample
  mostAppearingNames <- lapply(probability, base::names)

  # Store the one with the highest probability score into the bestFit dataframe
  for (j in seq(1:length(mostAppearingNames))) {
    numberPositions <- as.numeric(positions[j])
    probabilityScores[j] <- probability[[j]][numberPositions]
    bestFit[j,] <- mostAppearingNames[[j]][numberPositions]
  }


  # Store the bestFit, the Originalcall and the accompanying probability score within the final dataframe.

  if (addOriginalCall == T) {
    # Look at the original calls for each test sample
    originalCall <- metaDataRef[rownames(result),classColumn]

    classifications <- cbind(predict = bestFit,
                             originalCall = originalCall,
                             probability = probabilityScores)
  } else {
    classifications <- cbind(predict = bestFit,
                             probability = probabilityScores)
  }


  if (nrow(result) == 1) {
    classifications <- classifications[1, , drop = F]

  }
  # Make sure that the classifications still have their accompanying biomaterial_id
  rownames(classifications) <- rownames(result)

  classificationList <- list(classifications = classifications,
                             probabilityList = probability,
                             metaDataRef = metaDataRef
  )

  return(classificationList)

}
