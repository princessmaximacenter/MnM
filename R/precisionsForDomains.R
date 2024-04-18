#' Compare results between samples from different domains
#'
#' @param classColumn Column in the original metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Column in the original metadata file that contains the tumor type labels.
#' @param domainColumn Column in the original metadata file that contains the domain labels.
#' @param minorityDir Directory in which the minority model(s) are stored.
#' @param majorityDir Directory in which the majority model(s) are stored.
#' @param subtype Do you want to obtain the predictions on the tumor subtype classification level?
#' @param metaDataTest Metadata file containing the links between the patients and the tumor (sub)type diagnosis within the test set.
#' This file needs to have the same column names for the tumor subtype labels (classColumn),
#' tumor type labels (higherClassColumn) and domain labels (domainColumn).
#' @param probabilityThreshold What is the probability score threshold you would like to use to call a classification 'confident' for M&M?
#' @param metaDataRef Metadata file containing the links between the patients and the tumor (sub)type diagnosis.
#'
#' @return Dataframe containing the mean accuracy ($meanAccuracy),
#' precision ($meanPrecision), and recall ($meanRecall),
#' together with their standard deviations ($sdAccuracy,
#' $sdPrecision, $sdRecall), and the total number of samples associated to the domain group ($nSamples)
#' @export
#'
precisionsForDomains <- function(
          classColumn,
          higherClassColumn,
          domainColumn,
          minorityDir,
          majorityDir,
          subtype,
          metaDataTest,
          probabilityThreshold,
          metaDataRef) {

  `%notin%` <- Negate(`%in%`)
  if ("minorityClassifierResult.rds" %notin% list.files(minorityDir)) {
    crossValidation <- T
    allDirsMinority <- list.dirs(minorityDir, recursive = F)
    allDirsMajority <- list.dirs(majorityDir, recursive = F)
    selectedDirsMinority <- allDirsMinority[grep("seed", allDirsMinority)]
    selectedDirsMajority <- allDirsMajority[grep("seed", allDirsMajority)]

    print(paste0("Found ",length(selectedDirsMajority), " directories with different cross-validation runs.",
                 " Calculating average performance values for all combined."))
    if (length(selectedDirsMinority) != length(selectedDirsMajority)) {
      stop("The number of models for the minority and majority classifier are not the same.
         Please check your models within the minorityDir and majorityDir that the
         same seeds have been used for the generation of a minority and a majority classifier.")
    } else if (!all.equal(selectedDirsMajority, selectedDirsMinority) ) {
      stop("Please make sure you run the crossvalidation with the same seed for complementary classifications,
         and store them in the same directory.")
    }
  } else {
    selectedDirsMajority <- majorityDir
    crossValidation <- F
  }

  for (i in seq(1:length(selectedDirsMajority))) {
    if (crossValidation == T) {
      minorityDoc <- paste0(selectedDirsMinority[i],"/crossValidationMinorityResults.rds")
      majorityDoc <- paste0(selectedDirsMajority[i],"/crossValidationMajorityResults.rds")

    } else {
      minorityDoc <- paste0(minorityDir, "/minorityClassifierResult.rds")
      majorityDoc <- paste0(majorityDir, "/majorityClassifierResult.rds")
      metaData <- metaDataTest
      print("No original call found, adding it from metaDataTest")
    }
    minority <- readRDS(minorityDoc)
    majority <- readRDS(majorityDoc)

    metaData <- minority$metaDataRef
    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          subtype = subtype,
                                          classColumn = classColumn,
                                          higherClassColumn = higherClassColumn
    )

    predictionsMMFinal <- predictionsMMFinalList$predictionsMMFinal
    if ("originalCall" %in% colnames(predictionsMMFinal) ) {
      metaData <- minority$metaDataRef
      } else if (subtype == F) {
        predictionsMMFinal$originalCall <- metaData[rownames(predictionsMMFinal), higherClassColumn]
      } else {
        predictionsMMFinal$originalCall <- metaData[rownames(predictionsMMFinal), classColumn]
      }


    ourDomains <- unique(minority$metaDataRef[,domainColumn])


    for (j in seq(1:length(ourDomains))) {

      samplesDomain <- metaData %>% filter(Domain == ourDomains[j]) %>% rownames(.)
      nSamples <- length(samplesDomain)
      predictionsMMDomain <- predictionsMMFinal %>% filter(rownames(.) %in% samplesDomain,
                                                           probability1 > probabilityThreshold)
      correct <- predictionsMMDomain %>% filter(originalCall == predict) %>% nrow()
      correctAll <- predictionsMMFinal %>% filter(rownames(.) %in% samplesDomain,
                                                  originalCall == predict) %>% nrow()
      accuracy <-  correctAll / length(samplesDomain)
      precision <- correct / nrow(predictionsMMDomain)
      recall <- nrow(predictionsMMDomain) / length(samplesDomain)

      ourDF <- data.frame(Domain = ourDomains[j],
                          nSamples = nSamples,
                          accuracy = accuracy,
                          precision = precision,
                          recall = recall)
      if (j == 1) {
        totalDF <- ourDF
      } else {
        totalDF <- rbind(totalDF, ourDF)

      }
    }

      nSamples <- nrow(predictionsMMFinal)
      predictionsMMFinalFiltered <- predictionsMMFinal %>% filter( probability1 > probabilityThreshold)
      correct <- predictionsMMFinalFiltered %>% filter(originalCall == predict) %>% nrow()
      correctAll <- predictionsMMFinal %>% filter(originalCall == predict) %>% nrow()
      accuracy <- correctAll / nrow(predictionsMMFinal)
      precision <- correct / nrow(predictionsMMFinalFiltered)
      recall <- nrow(predictionsMMFinalFiltered) / nrow(predictionsMMFinal)

      ourDF <- data.frame(Domain = "All",
                          nSamples = nSamples,
                          accuracy = accuracy,
                          precision = precision,
                          recall = recall)

      totalDF <- rbind(totalDF, ourDF)




      totalDF$seed <- i
    if (i == 1) {
      accuracyDF <- totalDF
    } else {
      accuracyDF <- rbind(accuracyDF, totalDF)
    }
  }

  meanNumbers <- accuracyDF %>% group_by(Domain) %>%
    summarise(
      nSamples = mean(nSamples),
      meanAccuracy = mean(accuracy),
      meanPrecision = mean(precision),
      meanRecall = mean(recall),
      sdAccuracy = sd(accuracy),
      sdPrecision = sd(precision),
      sdRecall = sd(recall)
    )

 return(meanNumbers)

}
