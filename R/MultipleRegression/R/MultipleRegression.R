# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("p.value", "estimate", "2.5 %", "97.5 %"))

#' Function for running a multiple linear regression model with listed variables
#'
#' @param input_data "mids" object (imputed dataset)
#' @param mlr_var_list A character vector containing the names of exposures to include in the multiple model
#' @param variable_list A dataframe with short and full names of exposure variables
#' @param conf_list A character vector specifying names of model adjustment factors (confounders)
#' @param sdq A string defining which SDQ score will be fitted
#'
#' @importFrom MASS glm.nb
#' @import dplyr
#' @importFrom stringr str_c
#' @importFrom mice pool
#' @return A dataframe with regression results for all exposures in the multiple model
#' @export
#'
MultipleRegression <- function(input_data, mlr_var_list, variable_list, conf_list, sdq) {

  # Run multiple regression for sdq ~ exposure + co-exposures + covariates
  fit_mlr <- with(input_data, MASS::glm.nb(as.formula(paste(sdq, paste(append(mlr_var_list, conf_list), collapse = " + "), sep = " ~ "))))

  # pool regression results
  pool_fit <- fit_mlr %>%
    mice::pool() %>%
    summary(conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::rename(Variable_name = term) %>%
    dplyr::filter(Variable_name %in% mlr_var_list) %>%
    dplyr::mutate(p_value = format(round(p.value, 3), nsmall = 3),
                  estimate = format(round(estimate, 2), nsmall = 2),
                  conf_low = format(round(`2.5 %`, 2), nsmall = 2),
                  conf_high = format(round(`97.5 %`, 2), nsmall = 2),
                  estCI = stringr::str_c(estimate, " (", conf_low, "; ", conf_high, ")", sep = "")) %>%
    merge(variable_list, ., by = "Variable_name") %>%
    dplyr::select(Variable_name, estCI, p_value)

  return(pool_fit)
}
