# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("HelixID", "cohort", "sex", "h_sdq_external"))

#' Create a one (number 20) complete imputed dataset (from mids object)
#'
#' @param input_data A mids object
#' @param variable_list Exposure variables to be analysed
#' @param indx A scalar from 1 to n of imputed datasets, choosing which dataset will be returned
#'
#' @return A 'random' complete dataset (number 20)
#' @export
#' @import mice
#' @import dplyr
#' @importFrom tidyselect all_of

Mids2SingleComplete <- function(input_data, variable_list, indx = 20) {

  # Create one complete imputed dataset (for example basing on imputed dataset 20)
  all_imps <- mice::complete(input_data, indx) %>%

    # Select variables of interest plus covariates, SDQ cores and imp number
    dplyr::select(HelixID, cohort, h_sdq_external:sex, tidyselect::all_of(variable_list$Variable_name)) %>%
    droplevels()

  return(all_imps)
}
