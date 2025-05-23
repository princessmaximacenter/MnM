checkFormatInputData <- function(sampleColumn,
                                 classColumn,
                                 higherClassColumn,
                                 domainColumn,
                                 metaDataRef,
                                 countDataRef,
                                 outputDir = "NA",
                                 saveModel) {

  if (sampleColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the sample IDs is not present within metaDataRef. Please check the sampleColumn.")
  } else if (classColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the tumor subtype labels is not present within metaDataRef. Please check the classColumn")
  } else if (higherClassColumn %notin% base::colnames(metaDataRef)){
    base::stop("The column you specified for the tumor type labels is not present within metaDataRef. Please check the higherClassColumn")
  } else if (domainColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the tumor domain labels is not present within metaDataRef. Please check the domainColumn")
  }

  # Make sure the metadata and count data are in the right format and same order
  if (base::nrow(metaDataRef) != base::ncol(countDataRef)) {
    base::stop("The number of samples do not match between the metadata and the count data. Please make sure you include all same samples in both objects.")
  } else if (base::all(metaDataRef[, sampleColumn] %notin% base::colnames(countDataRef))) {
    base::stop("Your input data is not as required. Please make sure your sample IDs are stored in the sampleColumn, and in the column names of the count data")
  }


  if (base::is.numeric(countDataRef) != T) {
    base::stop("Your input data is not as required. Please make sure your countDataRef object only contains numerical count data and is a matrix.")

  }

  # Include a statement to store the classColumn, higherClassColumn and domainColumn
  base::cat(base::paste0("The column used for tumor subtypes labels within the metadata, used for model training purposes, is: ",
                           classColumn, '\nThis column contains values such as: \n'))
  base::cat(paste0(base::unique(metaDataRef[,classColumn])[1:3]), "\n")

  base::cat(base::paste0("\nThe column used for tumor type labels within the metadata, is: ",
                         higherClassColumn,'\nThis column contains values such as: \t'))
  base::cat(paste0(base::unique(metaDataRef[,higherClassColumn])[1:3]), "\n\n")

  base::cat(base::paste0("\nThe column used for tumor domain labels within the metadata, is: ",
                         domainColumn, '\nThis column contains values such as: \t'))
  base::cat(base::unique(metaDataRef[,domainColumn])[1:3])
  base::cat(paste0("\n\nIf any of these are incorrect, specify a different 'classColumn' (subtype),",
            "\n'higherClassColumn' (tumor type) or 'domainColumn' (domain) to function as labels.\n\n"))

  if (saveModel == T & !base::dir.exists(outputDir)) {
    checkDirectory <- base::tryCatch(base::dir.create(outputDir))
    if (checkDirectory == F) {
      base::stop(base::paste0("The directory you want the classification to be saved in cannot be created due to an error in the directory path.",
                              " Please check the spelling of your specified outputDir - it is probable the parent-directory does not exist."))
    }
  }

}


checkFormatTestData <- function(countDataNew,
                                countDataRef,
                                outputDir,
                                saveModel) {

  if (is.null(countDataRef)) {
    base::stop("You probably are using an old version of the model that is no longer compatible with MnM.\n\nPlease download the latest version from https://zenodo.org/records/14167359. ")
  }
  # Check whether the genes are within the rows of countDataNew
  geneOverlap <- sum(rownames(countDataNew) %in% rownames(countDataRef))

   if(geneOverlap == 0) {
   base::stop("Your input data is not as required. Please make sure your genes are in the rownames of countDataNew and are in the format of HGNC gene names.")
 } else if(geneOverlap < 0.6 * nrow(countDataRef)) {
  base::cat("Please note that there is less than 60% overlap between the genes within the reference cohort and countDataNew.\n")

 }

  # Check whether counts are supplied within the RNA-seq counts
  if (base::is.numeric(countDataNew) != T) {
    base::stop("Your input data is not as required. Please make sure your countDataNew object only contains numerical count data and is a matrix.")
  }

  # Generate the directory, if not possible abort
  if (saveModel == T & !base::dir.exists(outputDir)) {
    checkDirectory <- base::tryCatch(base::dir.create(outputDir))
    if (checkDirectory == F) {
      base::stop("The directory you want the classification to be saved in cannot be created due to an error in the directory path. Please check the spelling of your specified outputDir.")
    }
  }

}

