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
#' @param classColumn Column in the metadata file that contains the tumor (sub)type labels.
#' @param domainColumn Column in the metadata file that contains the tumor domain labels.
#' @param abbreviationTumorType Dataframe containing the links between the tumor (sub)type,
#' the abbreviation required in the plot, and the domain.
#' @param proteinFile In which directory can we find the file specifying the names of protein-coding genes within our dataset?
#' @param whichSeed For reproducibility, the seed can be specified with this parameter.
#'
#' @return List containing the UMAP-transformed datapoints ($dataUMAP), the
#' umap-transformation model to help with transforming potential new samples
#' in the same space ($umapTransform), and the ribodepletion correction model ($riboModelList).
#' @export
#' @import uwot dplyr magrittr

createUMAPcohort <- function(countDataRef,
                             metaDataRef,
                             classColumn,
                             domainColumn,
                             abbreviationTumorType,
                             proteinFile = "~/surfdrive/Shared/Kemmeren group/Research_Projects/RNA_classification_FW/data/input/20230320_proteinCodingGenes_gencode31.csv",
                             whichSeed = 1) {


  proteinCodingGenes <- read.table(proteinFile, sep = "\t") %>%
    select(x) %>% deframe
  set.seed(whichSeed)
  riboModelList <- riboCorrectCounts(data = countDataRef,
                                     proteinCodingGenes = proteinCodingGenes,
                                     outputDir = ".",
                                     saveRiboModels = F
  )
  countDataRef <- riboModelList$counts

  # Log-transform data
  dataLogRef <- log(countDataRef +1) %>% t() %>% as.data.frame()

  metaDataJoined <- left_join(metaDataRef, abbreviationTumorType[,c(classColumn, "abbreviation")])
  rownames(metaDataJoined) <- rownames(metaDataRef)
  metaDataRef <- metaDataJoined

  dataLogRef$abbreviation <- metaDataRef[rownames(dataLogRef), "abbreviation"]
  dataLogRef$subclass <- metaDataRef[rownames(dataLogRef), classColumn]
  dataLogRef$Domain <- metaDataRef[rownames(dataLogRef), domainColumn]


  set.seed(whichSeed)
  dataLogUMAP <- dataLogRef %>%
    select(where(is.numeric)) %>%
    uwot::umap(., ret_model = T)

  colnames(dataLogUMAP$embedding) <- c("UMAP1", "UMAP2")

  dataUMAP <- dataLogUMAP$embedding %>%
    as.data.frame()

  dataUMAP <- cbind(dataUMAP, dataLogRef[, c("subclass","Domain", "abbreviation"), drop = F])
  dataUMAPList <- list(dataUMAP = dataUMAP,
       umapTransform = dataLogUMAP,
       riboModelList = riboModelList)
  return(dataUMAPList)

}
