#' Function used to process RNA-Seq count data from a text file and save as a CSV
#' This data is obtained from a Diagnostics RNA-Seq based Pipeline
#'
#' This function processes RNA-Seq count data from a specified input file, summarizes the
#' counts by gene, and outputs the summarized data as a CSV file. The data is grouped by
#' gene name, and CPM values are calculated as the sum of counts per gene.
#'
#' @param input_file_path String representing the path to the input RNA-Seq count data file.
#' @param output_file_path String representing the path where the output CSV file will be saved.
#' @param sample_name String representing the name of the column for which CPM values are calculated, this is the sample_name

#' @return None. The function saves the processed data as a CSV file to the specified path.
#'
processCountData <- function(input_file_path, output_file_path, sample_name) {
  # Read the input file lines
  lines <- readLines(input_file_path)
  
  # Extract the header line that starts with '# ID'
  header_line <- grep("^# ID", lines, value = TRUE)
  
  # Split the header into column names
  col_names <- unlist(strsplit(sub("^# ", "", header_line), "\t"))
  
  # Read the data file, skipping lines starting with '#'
  data <- read.table(input_file_path, header = FALSE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
  
  # Assign column names to the data
  colnames(data) <- col_names
  
  # Group data by GeneName and calculate the sum of CPM values
  data_summarised <- data %>%
    dplyr::group_by(GeneName) %>%
    dplyr::summarise(CPM = sum(CPM, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  # Set GeneName as rownames
  data_summarised <- data_summarised %>%
    tibble::column_to_rownames("GeneName")
  
  # Rename the CPM column for clarity
  colnames(data_summarised) <- c(sample_name)
  
  # Write the summarized data to a CSV file
  readr::write_csv(data_summarised, output_file_path)
}
