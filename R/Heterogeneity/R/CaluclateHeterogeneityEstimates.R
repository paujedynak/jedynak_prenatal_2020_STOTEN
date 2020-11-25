# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}


# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("est",
                         "se",
                         "V1",
                         "V2",
                         "V3",
                         "V4",
                         "Description",
                         "Family",
                         "Exposure",
                         "Variable_name"))

#' Calculate heterogeneity estimates - CODE and COMMENTS come from Xavier Basagaña
#'
#' @param data A complete stacked dataset (mice object), including .imp = 0
#' @param outcome A string defining an outcome ("h_sdq_external" or "h_sdq_ternal")
#' @param conf_list A vector containing all confounder variable names
#' @param variable_list A data.frame containing a list of exposures for which regressions will be run, first column must contain exposure names
#' @param strat_variable A string identifying the name of the variable to stratify on
#'
#' @return A dataframe with heterogeneity estimates for all exposures of interest
#' @export
#' @import metaplus
#' @import dplyr
#' @importFrom tidyr spread

CaluclateHeterogeneityEstimates <- function(data, outcome, conf_list, variable_list, strat_variable) {

  selected_exposures <- variable_list$Variable_name
  estimates <- .CalculateI(data, outcome, conf_list, selected_exposures, strat_variable)

  est_result <- data.frame()

  for (k in seq_along(selected_exposures)) {

    temp <- estimates$RESv[[selected_exposures[k]]]
    mag_meta <- metaplus::metaplus(yi = est, sei = se, slab = rownames(temp), data = data.frame(temp), plotci = FALSE)
    mag_meta$ci_lb <- mag_meta$yi - (1.96 *  mag_meta$sei)
    mag_meta$ci_ub <- mag_meta$yi + (1.96 *  mag_meta$sei)

    mag_meta$est_ci <- paste0(format(round(exp(mag_meta$yi), 2), nsmall = 2),
                              " (",
                              format(round(exp(mag_meta$ci_lb), 2), nsmall = 2),
                              ";",
                              format(round(exp(mag_meta$ci_ub), 2), nsmall = 2),
                              ")")

    if (estimates$RES[[k]]$I2/100 < 0.001) {
      I <- paste0("<", formatC(0.001, digits = 4))

    } else {
      I <- format(round(estimates$RES[[k]]$I2/100, 3), nsmall = 3)
      I <- formatC(I, digits = 4)
    }

    mag_meta$I2 <- I

    estim <- as.data.frame(cbind(selected_exposures[k], mag_meta$slab, mag_meta$est_ci, mag_meta$I2))

    est_result <- rbind(est_result, estim)
  }

  est_result <- tidyr::spread(data = est_result, key = V2, value = V3) %>%
    dplyr::rename(Variable_name = V1, I2 = V4)

  est_result <- est_result[, c(1, 3:7, 2)] %>%
    merge(dplyr::select(variable_list, -Description), ., by = "Variable_name") %>%
    dplyr::arrange(Family, Exposure)%>%
    dplyr::select(-Variable_name, Family)

  return(est_result)
}
