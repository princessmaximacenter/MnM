#' UMAP-transform reference cohort
#'
#' Function to create the data required for plotting an UMAP of the reference cohort,
#' and paving the way for transforming new samples into the same sample space by
#' saving the transformation-model and ribo-depletion correction model.
#'
#' @param countDataRef Matrix containing the RNA-transcript per million data. Patients are in the columns,
#' different genes in the rows.
#' @param metaDataRef Metadata file containing the links between the patients and
#' the tumor (sub)type diagnosis within the reference cohort.
#' @param classColumn Column in the metadata file that contains the tumor subtype labels.
#' @param domainColumn Column in the metadata file that contains the tumor domain labels.
#' @param abbreviations Dataframe containing the links between the tumor (sub)type,
#' the abbreviation required in the plot, and the domain.
#' @param proteinFile In which directory can we find the file specifying the names of protein-coding genes within our dataset?
#' @param whichSeed For reproducibility, the seed can be specified with this parameter.
#' @param higherClassColumn Column in the metadata file that contains the tumor type labels.
#' @param correctRibo Do you want to perform a correction for the ribodepletion protocol on your dataset?
#' @param subtype Do you want to both get the tumor type and subtype labels within your UMAP object?
#' Default = FALSE, giving only the tumor type labels with the associated abbreviations.
#'
#' @return List containing the UMAP-transformed datapoints ($dataUMAP), and the ribodepletion correction model ($riboModelList).
#' @export
#' @import umap dplyr magrittr

createUMAPcohort <- function(countDataRef,
                             metaDataRef,
                             classColumn,
                             higherClassColumn,
                             domainColumn,
                             correctRibo = T,
                             abbreviations,
                             subtype =F,
                             proteinFile = "~/surfdrive/Shared/Kemmeren group/Research_Projects/RNA_classification_FW/data/input/20230320_proteinCodingGenes_gencode31.csv",
                             whichSeed = 1) {

  if (correctRibo == T) {
  proteinCodingGenes <- read.table(proteinFile, sep = "\t") %>%
    select(x) %>% deframe
  set.seed(whichSeed)
  riboModelList <- riboCorrectCounts(data = countDataRef,
                                     proteinCodingGenes = proteinCodingGenes,
                                     outputDir = ".",
                                     saveRiboModels = F
  )
  countDataRef <- riboModelList$counts
  }
  # Log-transform data
  dataLogRef <- log(countDataRef +1) %>% t() %>% as.data.frame()
  if (subtype == T) {
  abbreviations %<>% filter(!!sym(classColumn) %in% unique(metaDataRef[, classColumn]))

  metaDataJoined <- left_join(metaDataRef, abbreviations[,c(classColumn, "abbreviation")])
  } else {
    abbreviations %<>% filter(!!sym(higherClassColumn) %in% unique(metaDataRef[, higherClassColumn]))

    metaDataJoined <- left_join(metaDataRef, abbreviations[,c(higherClassColumn, "abbreviation")])

  }
  rownames(metaDataJoined) <- rownames(metaDataRef)
  metaDataRef <- metaDataJoined

  set.seed(whichSeed)
  dataLogUMAP <- dataLogRef %>%
    select(where(is.numeric)) %>%
    umap::umap()
  colnames(dataLogUMAP$layout) <- c("UMAP1", "UMAP2")

  dataUMAP <- dataLogUMAP$layout %>%
    as.data.frame()

  dataUMAP$abbreviation <- metaDataRef[rownames(dataUMAP), "abbreviation"]

  dataUMAP$Domain <- metaDataRef[rownames(dataUMAP), domainColumn]

if (subtype == T) {
  dataUMAP$subclass <- metaDataRef[rownames(dataUMAP), higherClassColumn]
  dataUMAP$subtype <- metaDataRef[rownames(dataUMAP), classColumn]

} else {
  dataUMAP$subclass <- metaDataRef[rownames(dataUMAP), classColumn]
}


  #dataUMAP <- cbind(dataUMAP, dataLogRef[, c("subclass","Domain", "abbreviation"), drop = F])
  dataUMAPList <- list(dataUMAP = dataUMAP)

  if (correctRibo == T) {
    dataUMAPList$riboModelList <- riboModelList
  }
  return(dataUMAPList)

}
