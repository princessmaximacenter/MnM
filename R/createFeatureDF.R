#' Title
#'
#' @param ml_models
#' @param df_of_interest
#' @param column
#' @param no_features
#'
#' @return
#' @export
#'
createFeatureDF <- function(ml_models, df_of_interest, column, no_features) {
  ##Create an empty dataframe
  empty_df <- data.frame()

  #Recursively extract column values from ModelList
  for (i in seq_along(ml_models)) {

    columnValues <- as.data.frame(ml_models[[i]][[df_of_interest]]) %>%
      select(all_of(column)) %>%
      arrange(desc(column)) %>%
      slice(1:no_features) %>%
      rownames_to_column(var = "Gene") %>%
      mutate(model = i) %>%
      rename(value = column)

    #Append to Dataframe
    empty_df <- rbind(empty_df, columnValues)

  }

  #Create create a wide DF for comparing features
  wide_df <- empty_df %>% pivot_wider(names_from = "model")

  #Store columnName as rownames
  wide_df <- wide_df %>% remove_rownames %>% column_to_rownames(var="Gene")

  return(wide_df)
}
