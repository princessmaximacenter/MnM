#' Title
#'
#' @param refCohort
#' @param inputDir
#'
#' @return
#' @export
#'
#' @examples
removeDuplicates <- function(refCohort, inputDir) {

  metaData <- refCohort$metaData
  countData <- refCohort$counts

  ###CSR Cohort
  CSR_biomaterials <- read.csv(paste0(inputDir, "DPIL_CSR_Biomaterial_RNA_20230302.csv"), sep = ";", header = T) %>% rename( "biosource_id" = "src_biosource_id") #%>% mutate_all(na_if,"")
  CSR_biosource <- read.csv(paste0(inputDir, "DPIL_CSR_Biosource_20230302.csv"), sep = ";", header = T) #%>% mutate_all(na_if,"")
  CSR_diagnosis <- read.csv(paste0(inputDir, "HiX_CSR_Diagnosis_20230302.csv"), sep = ";", header = T) #%>% mutate_all(na_if,"")
  CSR_individual <- read.csv(paste0(inputDir, "HiX_CSR_Individual_20230302.csv"), sep = ";", header = T) #%>% mutate_all(na_if,"")

  CSR_metaData <- CSR_biomaterials %>% left_join(CSR_biosource %>% select(individual_id, biosource_id, diagnosis_id, disease_status, biosource_date, disease_status), by = c("individual_id", "biosource_id")) %>%
    left_join(CSR_diagnosis %>% select(individual_id, diagnosis_id, diagnosis_date, tumor_type_label), by = c("individual_id", "diagnosis_id")) %>%
    left_join(CSR_individual %>% select(individual_id, IC_given_date, IC_withdrawn_date), by = "individual_id") %>%
    unique %>% mutate_all(na_if,"")


  ###Removing Duplicates
  PMCIDs <- metaData[duplicated(metaData$PMCID),"PMCID"]
  duplicateEntries <- metaData %>% filter(PMCID %in% PMCIDs)
  duplicateEntries$Patient <- rownames(duplicateEntries)

  duplicateEntries <- left_join(duplicateEntries, CSR_metaData[, c("biosource_date", "biomaterial_id")],
                                by = c("Patient" = "biomaterial_id"))

  rownames(duplicateEntries) <- duplicateEntries$Patient

  #rm(nestedDF)
  #1 sample per patient
  nestedDF <- duplicateEntries %>% group_by(PMCID, Disease_sub_class) %>%
    mutate(biosource_date = as.Date(biosource_date, "%d/%m/%Y")) %>% nest

  #group first by patient and then disease_subclass (we want to consider mismatches)
  #Look to trecode data to fill in missing values for biosource_dates

  keepSampleDuplicates <- c()

  for (i in seq(1:length(nestedDF$data))) {
    #biosource_date < biomaterial_date < diagnosis_date -> for arranging
    nestedDF$data[[i]] %<>% arrange(biosource_date) #change to biosource_date
    keepSample <- nestedDF$data[[i]][1,"Patient", drop = T]
    keepSampleDuplicates <- c(keepSampleDuplicates, keepSample)
  }
  removeEntries <- duplicateEntries %>% filter(rownames(.) %notin% keepSampleDuplicates)

  newMetaData <- metaData %>% filter(rownames(.) %notin% rownames(removeEntries))

  all_RNA_seq_biomaterial_IDs <- countData %>%
    colnames()

  selected_RNA_seq_biomaterial_IDs <- rownames(newMetaData)

  newCountData <- countData[, all_RNA_seq_biomaterial_IDs %in% selected_RNA_seq_biomaterial_IDs]

  refCohort$counts <- newCountData
  refCohort$metaData <- newMetaData

  return(refCohort)
}
