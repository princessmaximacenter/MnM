#' Combine individual model results into final classification
#'
#' Convert individual classifications to an object containing the final per-sample classifications,
#' per-label probabilities and metadata that the classifications were based on.
#'
#' @param result Dataframe containing the classification labels from each model for each sample.
#' @param metaDataRef Metadata file containing the links between the samples and the tumor domain, type and subtype diagnosis.
#' @param addOriginalCall Would you like to add the original tumor subtype label from the metadata (TRUE)?
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#'
#' @return List containing the final per-sample classifications ($classifications),
#' per-label probabilities ($probabilityList) and metadata that the classifications were based on ($metaDataRef).

convertResultToClassification <- function(result,
         metaDataRef,
         addOriginalCall,
         classColumn = NA) {

  if (base::nrow(result) == 1 || typeof(base::apply(result, 1, base::table)) == 'integer') {
    randomVector <- base::paste0("fake", 1:base::ncol(result)) %>%
      base::as.data.frame() %>%
      base::t() %>% base::as.data.frame()
    base::colnames(randomVector) <- base::colnames(result)
    result1 <- base::rbind(result, randomVector)
    probability <-  base::apply(result1, 1, base::table)
    probability <- probability[-base::length(probability)]
  } else {
    # Find out how often a certain tumor type prediction is made for a specific sample
    probability <- base::apply(result, 1, base::table)
  }
  # Locate the position of the highest probability
  positions <- base::lapply(probability, base::which.max)
  positions <- base::unlist(positions)

  bestFit <- base::data.frame(predict = base::rep(NA, times = base::length(result$fold1)))
  probabilityScores <- base::vector()

  # Extract the different calls being made for each sample
  mostAppearingNames <- base::lapply(probability, base::names)

  # Store the one with the highest probability score into the bestFit dataframe
  for (j in base::seq(1:base::length(mostAppearingNames))) {
    numberPositions <- base::as.numeric(positions[j])
    probabilityScores[j] <- probability[[j]][numberPositions]
    bestFit[j,] <- mostAppearingNames[[j]][numberPositions]
  }


  # Store the bestFit, the Originalcall and the accompanying probability score within the final dataframe.

  if (addOriginalCall == T) {
    # Look at the original calls for each test sample
    originalCall <- metaDataRef[base::rownames(result),classColumn]

    classifications <- base::cbind(predict = bestFit,
                             originalCall = originalCall,
                             probability = probabilityScores)
  } else {
    classifications <- base::cbind(predict = bestFit,
                             probability = probabilityScores)
  }


  if (base::nrow(result) == 1) {
    classifications <- classifications[1, , drop = F]

  }
  # Make sure that the classifications still have their accompanying biomaterial_id
  base::rownames(classifications) <- base::rownames(result)

  classificationList <- base::list(classifications = classifications,
                             probabilityList = probability,
                             metaDataRef = metaDataRef
  )

  return(classificationList)

}
