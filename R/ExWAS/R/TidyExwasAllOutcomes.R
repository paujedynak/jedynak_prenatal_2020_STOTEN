# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("2.5 %", "97.5 %", "estimate", "CI", "p_value", "p.value", "std.error", "statistic", "df", "Family",
                         "Exposure", "Estimate", "conf_low", "conf_high", "estCI", "p_value_FDR", "sub_cohort"))

#' Tidy the outcomes of regressions run for all exposures in a mids object
#'
#' @param input_data mids object (imputed dataset)
#' @param variable_list A data.frame containing a list of exposures for wich regressions will be run, first column must contain exposure names
#' @param exwas_result Regression output
#' @param subcohort A one element character vector indicating a stratification factor
#' @param save A character vector specifying if the files are going to be saved (default = TRUE, save files)
#' @param file_name File name (to which .csv or _plot.csv will be added)
#' @param path A string defining path for saving
#'
#' @return A data.frame with tidy regression result, ready for plotting
#' @import dplyr
#' @import utils
#' @importFrom here here
#' @importFrom mice complete
#' @importFrom stringr str_c
#' @importFrom stats cor
#'

.TidyExwasAllOutcomes <- function(input_data,
                                  variable_list,
                                  exwas_result,
                                  subcohort,
                                  save = FALSE,
                                  path = NULL,
                                  file_name = NULL) {

  # # Transform the results so they can be used for plotting and
  # # exported as .csv table

  # calculate effective number of
  # tests (Meff), to calculate FDR corrected p value)
  if (class(input_data) == "mids") {
    cormat <- mice::complete(data = input_data) %>%
      dplyr::select(tidyselect::all_of(variable_list$Variable_name)) %>%
      stats::cor(use = "pairwise.complete.obs")

  } else {

    # calculate correlation matrix
    cormat <- dplyr::select(input_data,
                            tidyselect::ends_with("_log2_iqr"),
                            tidyselect::ends_with("_ln_iqr")) %>%
      stats::cor(use = "pairwise.complete.obs")
  }

  # meff calculation based on Li 2012 method
  lambdas <- eigen(cormat)$values
  Meff <- ncol(cormat) - sum((lambdas > 1) * (lambdas - 1))  #30.06931 for 47 exposures, 32.69156 for 52 exposures

  # Meff <- poolR::meff(cormat, method = "liji") # 35 exposures

  # Add p value corrected for multiple comparison for
  # correlated explanatory variables, Li 2012 method
  exwas_result <- exwas_result %>%
    dplyr::mutate(CI = stringr::str_c(format(round(`2.5 %`, 2), nsmall = 2),
                                      format(round(`97.5 %`, 2), nsmall = 2),
                                      sep = "; "),
                  estCI = stringr::str_c(format(round(estimate, 2), nsmall = 2),
                                         " (", CI, ")",
                                         sep = ""),
                  p_value = format(round(p.value, 3), nsmall = 3),
                  p_value_FDR = format(round(pmin(1, Meff * p.value), 3), nsmall = 3)) %>% # Bonferroni correction
    dplyr::select(-c(std.error, statistic, df, p.value)) %>%
    dplyr::rename(Estimate = estimate, conf_low = "2.5 %", conf_high = "97.5 %")

  if (is.null(subcohort) == TRUE) {

    # Add exposure names
    tidy_exwas_res <- exwas_result %>%
      dplyr::mutate(Exposure = variable_list$Exposure,
                    Variable_name = term) %>%

      # add variable descriptions
      dplyr::left_join(dplyr::select(variable_list, Exposure, Family), ., by = "Exposure") %>%

      # Order the table by Exposure
      dplyr::arrange(Family, Exposure)

    # remove unnecessary variables
    pooled_fits_plot <- dplyr::select(tidy_exwas_res, Family, Exposure, Estimate, conf_low, conf_high, p_value)

    exwas_sel <- dplyr::select(tidy_exwas_res, Variable_name, Family, Exposure, estCI, p_value, p_value_FDR)

  } else {

    # Calculate for how many subcohorts the analysis was performed
    lev_num <- dplyr::select(exwas_result, sub_cohort) %>%
      unique() %>%
      nrow()

    # add variable descriptions
    tidy_exwas_res <- exwas_result %>%

      # Add a column with a name of subcohort for which the analysis was performed
      dplyr::mutate(Exposure = rep(variable_list$Exposure, times = lev_num)) %>%

      #Merge with exposures descriptions
      dplyr::left_join(dplyr::select(variable_list, Exposure, Family), ., by = "Exposure") %>%

      # Order the table by subcohort and Exposure
      dplyr::arrange(sub_cohort, Family, Exposure)

    # pooled_fits_plot <- dplyr::select(tidy_exwas_res, sub_cohort, Family, Exposure, Estimate, conf_low, conf_high, p_value)

    # Order the table by subcohort and Exposure
    exwas_sel <- dplyr::select(tidy_exwas_res, Variable_name, sub_cohort, Family, Exposure, estCI, p_value, p_value_FDR)
  }

  if (save == TRUE) {
    # export table result to a .csv file
    # utils::write.csv(pooled_fits_plot, file = here::here(path, paste0(file_name, "_plot.csv")))

    # export table result to a .csv file
    utils::write.csv(exwas_sel, file = here::here(path, paste0(file_name, ".csv")))
  }

  # return table ready to plot
  return(exwas_sel)
}
