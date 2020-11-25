# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("Variable_name", "expo_sdq_ext", "expo_sdq_int"))

#' Merges selected columns from 2 datasets
#'
#' @param input_data_ext A dataframe with SDQ externalising ExWAS results
#' @param input_data_int A dataframe with SDQ internalising ExWAS results
#' @param expo_sdq_ext A character vector listing significant hits in ExWAS for SDQ external
#' @param expo_sdq_int A character vector listing significant hits in ExWAS for SDQ internal
#'
#' @importFrom dplyr filter
#' @importFrom tidyselect all_of
#' @return A dataframe consisting of 2 merged datasets
#' @export

MergeDatasets <- function(input_data_ext, input_data_int, expo_sdq_ext, expo_sdq_int) {

dataset <- dplyr::filter(input_data_ext, Variable_name %in% tidyselect::all_of(expo_sdq_ext)) %>%
  rbind(dplyr::filter(input_data_int, Variable_name %in% tidyselect::all_of(expo_sdq_int)))

return(dataset)
}
