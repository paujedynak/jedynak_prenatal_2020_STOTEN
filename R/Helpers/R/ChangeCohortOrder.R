#' Change order of the cohorts in a dataset to "BiB", "EDEN", "INMA", "KANC", "RHEA"
#'
#' @param input_data A dataframe
#'
#' @importFrom dplyr mutate
#' @return A dataframe with new order of the cohorts
#' @export
#'
ChangeCohortOrder <- function(input_data) {

  output_data <- input_data %>%
    dplyr::mutate(cohort = factor(cohort,
                                  levels = sort(levels(cohort)),
                                  labels = c("BiB", "EDEN", "INMA", "KANC", "RHEA")))

  return(output_data)
}
