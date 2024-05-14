#' Compare results between samples from different domains
#'
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
          minorityDir,
          majorityDir,
          subtype,
          metaDataTest = NA,
          probabilityThreshold,
          metaDataRef) {

  `%notin%` <- base::Negate(`%in%`)
  if ("minorityClassifierResult.rds" %notin% base::list.files(minorityDir)) {
    crossValidation <- T
    allDirsMinority <- base::list.dirs(minorityDir, recursive = F)
    allDirsMajority <- base::list.dirs(majorityDir, recursive = F)
    selectedDirsMinority <- allDirsMinority[base::grep("seed", allDirsMinority)]
    selectedDirsMajority <- allDirsMajority[base::grep("seed", allDirsMajority)]

    base::print(base::paste0("Found ",base::length(selectedDirsMajority), " directories with different cross-validation runs.",
                 " Calculating average performance values for all combined."))
    if (base::length(selectedDirsMinority) != base::length(selectedDirsMajority)) {
      stop("The number of models for the minority and majority classifier are not the same.
         Please check your models within the minorityDir and majorityDir that the
         same seeds have been used for the generation of a minority and a majority classifier.")
    } else if (!base::all.equal(selectedDirsMajority, selectedDirsMinority) ) {
      stop("Please make sure you run the crossvalidation with the same seed for complementary classifications,
         and store them in the same directory.")
    }
  } else {
    selectedDirsMajority <- majorityDir
    crossValidation <- F
  }

  for (i in base::seq(1:base::length(selectedDirsMajority))) {
    if (crossValidation == T) {
      minorityDoc <- base::paste0(selectedDirsMinority[i],"/crossValidationMinorityResults.rds")
      majorityDoc <- base::paste0(selectedDirsMajority[i],"/crossValidationMajorityResults.rds")

    } else {
      minorityDoc <- base::paste0(minorityDir, "/minorityClassifierResult.rds")
      majorityDoc <- base::paste0(majorityDir, "/majorityClassifierResult.rds")
      metaData <- metaDataTest
    }
    minority <- base::readRDS(minorityDoc)
    majority <- base::readRDS(majorityDoc)

    classColumn <- minority$metaDataRun$classColumn
    higherClassColumn <- minority$metaDataRun$higherClassColumn
    domainColumn <- minority$metaDataRun$domainColumn

    if (!is.na(metaDataTest)[1]) {
      base::print("Checking the performance for the test set based on values provided in dataframe 'metaDataTest'.")

      if (classColumn %notin% colnames(metaDataTest)) {
        base::print("Please note that the wanted column for the tumor subtype labels cannot be found within 'metaDataTest'.")
        base::print(base::paste0("Either change the column with the tumor subtype labels to the name: ", classColumn))
        base::stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
      } else if (higherClassColumn %notin% base::colnames(metaDataTest)) {

        base::print("Please note that the wanted column for the tumor type labels cannot be found within 'metaDataTest'.")
        base::print(base::paste0("Either change the column with the tumor type labels to the name: ", higherClassColumn))
        base::stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")

      } else if (domainColumn  %notin% colnames(metaDataTest)) {
        base::print("Please note that the wanted column for the tumor domain labels cannot be found within 'metaDataTest'.")
        base::print(base::paste0("Either change the column with the tumor domain labels to the name: ", domainColumn))
        base::stop("Alternatively, use the function 'classColumns()' to substitute the class-column names within M&M's reference metadata to the name of your liking that's present within your own metaDataTest.")
      } else {
        base::print(base::paste0("Found columns ", classColumn, ", ", higherClassColumn, ", and ", domainColumn, " within metaDataTest specifying the tumor subtype, type and domain."))

        base::print("No original call found, adding it from metaDataTest")
      }

    }
    predictionsMMFinalList <- integrateMM(minority = minority,
                                          majority = majority,
                                          subtype = subtype
    )

    predictionsMMFinal <- predictionsMMFinalList$predictionsMMFinal
    if ("originalCall" %in% base::colnames(predictionsMMFinal) ) {
      metaData <- minority$metaDataRef
      } else if (subtype == F) {
        predictionsMMFinal$originalCall <- metaData[base::rownames(predictionsMMFinal), higherClassColumn]
      } else {
        predictionsMMFinal$originalCall <- metaData[base::rownames(predictionsMMFinal), classColumn]
      }


    ourDomains <- base::unique(minority$metaDataRef[,domainColumn])


    for (j in base::seq(1:base::length(ourDomains))) {

      samplesDomain <- metaData %>% dplyr::filter(Domain == ourDomains[j]) %>% base::rownames(.)
      nSamples <- base::length(samplesDomain)
      predictionsMMDomain <- predictionsMMFinal %>% dplyr::filter(rownames(.) %in% samplesDomain,
                                                           probability1 > probabilityThreshold)
      correct <- predictionsMMDomain %>% dplyr::filter(originalCall == predict) %>% base::nrow()
      correctAll <- predictionsMMFinal %>% dplyr::filter(rownames(.) %in% samplesDomain,
                                                  originalCall == predict) %>% base::nrow()
      accuracy <-  correctAll / base::length(samplesDomain)
      precision <- correct / base::nrow(predictionsMMDomain)
      recall <- base::nrow(predictionsMMDomain) / base::length(samplesDomain)

      ourDF <- base::data.frame(Domain = ourDomains[j],
                          nSamples = nSamples,
                          accuracy = accuracy,
                          precision = precision,
                          recall = recall)
      if (j == 1) {
        totalDF <- ourDF
      } else {
        totalDF <- base::rbind(totalDF, ourDF)

      }
    }

      nSamples <- base::nrow(predictionsMMFinal)
      predictionsMMFinalFiltered <- predictionsMMFinal %>% dplyr::filter( probability1 > probabilityThreshold)
      correct <- predictionsMMFinalFiltered %>% dplyr::filter(originalCall == predict) %>% base::nrow()
      correctAll <- predictionsMMFinal %>% dplyr::filter(originalCall == predict) %>% base::nrow()
      accuracy <- correctAll / base::nrow(predictionsMMFinal)
      precision <- correct / base::nrow(predictionsMMFinalFiltered)
      recall <- base::nrow(predictionsMMFinalFiltered) / base::nrow(predictionsMMFinal)

      ourDF <- data.frame(Domain = "All",
                          nSamples = nSamples,
                          accuracy = accuracy,
                          precision = precision,
                          recall = recall)

      totalDF <- base::rbind(totalDF, ourDF)




      totalDF$seed <- i
    if (i == 1) {
      accuracyDF <- totalDF
    } else {
      accuracyDF <- rbind(accuracyDF, totalDF)
    }
  }

  meanNumbers <- accuracyDF %>% dplyr::group_by(Domain) %>%
    dplyr::summarise(
      nSamples = base::mean(nSamples),
      meanAccuracy = base::mean(accuracy),
      meanPrecision = base::mean(precision),
      meanRecall = base::mean(recall),
      sdAccuracy = stats::sd(accuracy),
      sdPrecision = stats::sd(precision),
      sdRecall = stats::sd(recall)
    )

 return(meanNumbers)

}
