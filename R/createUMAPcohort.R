#' UMAP-transform reference cohort
#'
#' Function to create the data required for plotting an UMAP of the reference cohort
#' and ribo-depletion correction model.
#'
#' Please note that the tumor domain, type and subtype all need to be specified within the metadata.
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data. Samples are in the columns,
#' different genes in the rows.
#' @param metaDataRef Metadata file containing the links between the samples and
#' the tumor (sub)type diagnosis within the reference cohort.
#' @param classColumn Name of column in the metadata file that contains the tumor subtype labels.
#' @param higherClassColumn Name of column in the metadata file that contains the tumor type labels.
#' @param domainColumn Name of column in the metadata file that contains the tumor domain labels.
#' @param sampleColumn Name of column in the metadata file that contains the sample identifiers.
#' @param abbreviations Optional. Dataframe containing the links between the tumor (sub)type,
#' the abbreviation required in the plot, and the domain.
#' @param proteinCodingGenes What are the names of the RNA-transcripts that stand for protein-coding genes within our dataset?
#' Please supply it as a vector. This is needed for ribo-depletion correction model.
#' @param whichSeed For reproducibility, the seed can be specified with this parameter.
#' @param correctRibo Do you want to perform a correction for the ribodepletion protocol on your dataset? Default is TRUE.
#' @param noAbbreviations Do you want to use abbreviations within the UMAP? If not, specify noAbbreviations = F. Default = TRUE.
#'
#' @return List containing the UMAP-transformed datapoints ($dataUMAP), the abbreviations used to name the UMAP-transformed datapoints, and the ribodepletion correction model ($riboModelList).
#' @export
#' @import umap

createUMAPcohort <- function(countDataRef,
                             metaDataRef,
                             classColumn,
                             higherClassColumn,
                             domainColumn,
                             sampleColumn,
                             correctRibo = F,
                             abbreviations = NA,
                             noAbbreviations = F,
                             proteinCodingGenes,
                             whichSeed = 1) {

  `%notin%` <- base::Negate(`%in%`)

  if (sampleColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the sample IDs is not present within metaDataRef. Please check the sampleColumn.")
  }  else if (classColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the tumor subtype labels is not present within metaDataRef. Please check the classColumn")
  } else if (higherClassColumn %notin% base::colnames(metaDataRef)){
    base::stop("The column you specified for the tumor type labels is not present within metaDataRef. Please check the higherClassColumn")
  } else if (domainColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the tumor domain labels is not present within metaDataRef. Please check the domainColumn")
  }

  if (base::nrow(metaDataRef) != base::ncol(countDataRef)) {
    base::stop("The number of samples do not match between the metadata and the count data. Please make sure you include all same samples in both objects.")
  }
  base::rownames(metaDataRef) <- metaDataRef[, sampleColumn]
  if (base::all(base::rownames(metaDataRef) %notin% base::colnames(countDataRef))) {
    base::stop("Your input data is not as required. Please make sure your sample IDs are within the row names of the metadata, and in the column names of the count data")
  }

  if (base::is.numeric(countDataRef) != T) {
    base::stop("Your input data is not as required. Please make sure your countDataRef object only contains numerical count data and is a matrix. Non-available measurements are not allowed.")

  }

  checkFormatInputData(sampleColumn = sampleColumn,
                       classColumn = classColumn,
                       higherClassColumn = higherClassColumn,
                       domainColumn = domainColumn,
                       metaDataRef = metaDataRef,
                       countDataRef = countDataRef,
                       saveModel = F
  )

  if (base::is.na(abbreviations)[1]) {

    abbreviations <- metaDataRef[, c(domainColumn, classColumn, higherClassColumn)] %>% base::unique()
    abbreviations$abbreviationSubtype <- abbreviations[,classColumn]
    abbreviations$abbreviationTumorType <- abbreviations[,higherClassColumn]

    base::print(base::paste("You have not supplied any abbreviations (abbreviations). If you would like to use abbreviations,",
                " please generate a dataframe with the following columns and abbreviations within abbreviationSubtype and abbreviationTumorType:"))
    base::print(abbreviations[1:4,])
  }

  if (correctRibo == T) {

  base::set.seed(whichSeed)
  riboModelList <- riboCorrectCounts(data = countDataRef,
                                     proteinCodingGenes = proteinCodingGenes,
                                     outputDir = ".",
                                     saveRiboModels = F
  )
  countDataRef <- riboModelList$counts
  }
  # Log-transform data
  dataLogRef <- base::log(countDataRef +1) %>% base::t() %>% base::as.data.frame()

  if (noAbbreviations == F) {
    abbreviations %<>% dplyr::filter(!!dplyr::sym(classColumn) %in% base::unique(metaDataRef[, classColumn]),
                                     !!dplyr::sym(higherClassColumn) %in% base::unique(metaDataRef[, higherClassColumn])
    )
    metaDataJoined <- dplyr::left_join(metaDataRef, abbreviations[,c(classColumn, "abbreviationTumorType", "abbreviationSubtype")])

    base::rownames(metaDataJoined) <- base::rownames(metaDataRef)
    metaDataRef <- metaDataJoined
  } else {
    metaDataRef$abbreviationTumorType <- metaDataRef[, higherClassColumn]
    metaDataRef$abbreviationSubtype <- metaDataRef[, classColumn]

  }


  base::set.seed(whichSeed)
  dataLogUMAP <- dataLogRef %>%
    dplyr::select(dplyr::where(is.numeric)) %>%
    umap::umap()
  base::colnames(dataLogUMAP$layout) <- c("UMAP1", "UMAP2")

  dataUMAP <- dataLogUMAP$layout %>%
    base::as.data.frame()

  if (noAbbreviations == F) {
    dataUMAP$abbreviationTumorType <- metaDataRef[base::rownames(dataUMAP), "abbreviationTumorType"]
    dataUMAP$abbreviationSubtype <- metaDataRef[base::rownames(dataUMAP), "abbreviationSubtype"]
  } else if (noAbbreviations == T) {
    dataUMAP$abbreviationTumorType <- metaDataRef[base::rownames(dataUMAP),higherClassColumn]
    dataUMAP$abbreviationSubtype <- metaDataRef[base::rownames(dataUMAP),classColumn]

  }

  dataUMAP$Domain <- metaDataRef[base::rownames(dataUMAP), domainColumn]


  dataUMAP$subclass <- metaDataRef[base::rownames(dataUMAP), higherClassColumn]
  dataUMAP$subtype <- metaDataRef[base::rownames(dataUMAP), classColumn]


  #dataUMAP <- cbind(dataUMAP, dataLogRef[, c("subclass","Domain", "abbreviation"), drop = F])
  dataUMAPList <- base::list(dataUMAP = dataUMAP,
                             abbreviations = abbreviations)

  if (correctRibo == T) {
    dataUMAPList$riboModelList <- riboModelList
  }
  return(dataUMAPList)

}
