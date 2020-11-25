# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("conf.low", "conf.high", "n_obs", "smoking", "h_sdq_external", "h_sdq_internal", "confounder", "coadj", "pvalue"))

#' Run an ExWAS on a dataset that does not need pooling afterwards (complete case dataset)
#'
#' @param input_data complete "mids" long format (stacked) imputed dataset or any dataset with no pooling required (complete cases for example)
#' @param variable_list A data.frame containing a list of exposures for wich regressions will be run, first column must contain exposure names
#' @param sdq A character vector indicating in the outcome SDQ score is external or internal
#' @param subcohort A character vector specifying stratification to be made (one cohort left)
#' @param save A character vector specifying if the files are going to be saved (default = TRUE, save files)
#' @param path Path for saving (default = NULL)
#' @param file_name File name (to which .csv or _plot.csv will be added, default = NULL)
#' @param conf_list A character vector defining a lsit of confounders used in the model
#'
#' @return Regression coefficients (IRR) with '95%' confidence intervals for complete dataset
#' @export
#' @import dplyr
#' @import utils
#' @importFrom broom tidy
#' @importFrom stats as.formula nobs
#' @importFrom here here
#' @importFrom MASS glm.nb


RunExwasCompleteCase <- function(input_data, variable_list, sdq = c("h_sdq_external", "h_sdq_internal"), conf_list, subcohort = NULL, save = FALSE, path = NULL, file_name = NULL) {

  exwas_result <- data.frame()

  # select one exposure
  for (item in variable_list$Variable_name) {

    # Create a formula for regression
    formula <- stats::as.formula(paste(sdq, "~", item, "+", paste(conf_list, collapse = " + ")))

    # run simple regression for sdq ~ exposure + covariates
    fit <- MASS::glm.nb(formula = formula, data = input_data)

    # tidy the result
    tidy_fit <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%

      # pick the row for an exposure only, ignore coef. for intercept and confounders
      .[2, ] %>%

      # add number of observations for the model
      dplyr::mutate(n_obs = stats::nobs(fit))

    # append the dataframe with results
    exwas_result <- rbind(exwas_result, tidy_fit)
  }

  clean_exwas_result <- exwas_result %>%
    dplyr::mutate(CI = stringr::str_c(format(round(conf.low, 2), nsmall = 2), format(round(conf.high, 2), nsmall = 2), sep = "; "),
                  estCI = stringr::str_c(format(round(estimate, 2), nsmall = 2), " (", CI, ")", sep = ""),
                  p_value = format(round(p.value, 3), nsmall = 3)) %>%
    dplyr::select(-c(std.error, statistic, p.value)) %>%
    dplyr::rename(Estimate = estimate, conf_low = conf.low, conf_high = conf.high, Variable_name = term)

  tidy_exwas_res <- clean_exwas_result %>%
    dplyr::mutate(Exposure = variable_list$Exposure) %>%

    # add variable descriptions
    dplyr::left_join(dplyr::select(variable_list, Exposure, Family), ., by = "Exposure") %>%

    # Order the table by Exposure
    dplyr::arrange(Family, Exposure)

  fits_plot <- dplyr::select(tidy_exwas_res, Family, Exposure, Estimate, conf_low, conf_high, p_value)

  # Order the table by subcohort and Exposure
  exwas_sel <- dplyr::select(tidy_exwas_res, Variable_name, estCI, p_value, n_obs)

  if (save == TRUE) {
    # export table result to a .csv file
    utils::write.csv(fits_plot, file = here::here(path, paste0(file_name, "_plot.csv")))

    # export table result to a .csv file
    utils::write.csv(exwas_sel, file = here::here(path, paste0(file_name, ".csv")))
  }

  return(exwas_sel)
}
