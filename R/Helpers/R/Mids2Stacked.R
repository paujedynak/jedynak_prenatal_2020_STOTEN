# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c(".imp", ".id", "cohort", "sex"))

#' Change mids object to a stacked data.frame
#'
#' @param input_data A mids object (imputed dataset)
#' @param variable_list A data.frame containing a list of exposures for which regressions will be run, first column must contain exposure names
#' @param include A logical to indicate whether the original data with the missing values should be included.
#'
#' @return A stacked imputed dataset (N x number of imputations)
#' @export
#' @import mice
#' @import dplyr
#' @importFrom tidyselect all_of


Mids2Stacked <- function(input_data, variable_list, include = FALSE) {

  # Create a complete, stacked dataset out of the 'mids' imputed dataset
  all_imps <- mice::complete(input_data, "long", include) %>%

    # Select variables of interest plus covariates, SDQ cores and imp number
    dplyr::select(.imp, .id, cohort:sex, tidyselect::all_of(variable_list$Variable_name)) %>%
    droplevels()

  return(all_imps)
}
