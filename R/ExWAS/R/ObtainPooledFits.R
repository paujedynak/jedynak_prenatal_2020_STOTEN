# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(".imp")


#' Pool the fits for imputed dataset regressions
#'
#' @param variable_list A data.frame containing a list of exposures for wich regressions will be run, first column must contain exposure names
#' @param conf_list A character vector specifying confounders to be added in the regression model
#' @param dataset A data.frame containing all stacked imputed datasets
#' @param sdq A character vector defining SDQ score
#'
#' @return A data.frame containing an untidy regression fits for all exposures
#' @import dplyr
#' @import utils
#' @importFrom mice pool
#' @importFrom stats as.formula
#' @importFrom MASS glm.nb
#'

.ObtainPooledFits <- function(dataset,
                              variable_list, sdq = c("h_sdq_external", "h_sdq_internal"),
                              conf_list) {

  i <- 1

  tidy_res <- data.frame()

  # Select one exposure at a time
  for (item in variable_list$Variable_name) {

    # Print progress
    print(paste0(i, "/", length(variable_list$Variable_name)))

    # Create a formula for regression
    formula <- stats::as.formula(paste(sdq, "~", item, "+", paste(conf_list, collapse = " + ")))

    # Run regression
    result <- dataset %>%

      # perform regression for each imputed dataset (group by
      # imputation number)
      dplyr::group_by(.imp) %>%

      # run simple regression for sdq ~ exposure + covariates
      dplyr::do(model = MASS::glm.nb(formula = formula, data = .)) %>%

      # transform the regression output to a list
      as.list() %>%

      # pick the list item containing numerical regression output
      # (coefficient and other statistics for the intercept,
      # confounders and exposure), ignore other results provided by
      # the glm.nb function
      .[[-1]] %>%

      # pool the regression outcomes for each imputed dataset
      mice::pool() %>%

      # extract the final coefficient with 95% confidence intervals
      summary(conf.int = TRUE, exponentiate = TRUE)


    # Tidy the result for one exposure
    pool_fit <- .TidyExwasSingleResult(result)

    # Merge the results for all the exposures
    tidy_res <- rbind(tidy_res, pool_fit)

    i <- i + 1
  }

  return(tidy_res)
}
