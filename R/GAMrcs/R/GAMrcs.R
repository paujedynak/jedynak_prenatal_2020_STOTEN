# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("Exposure", "Variable_name"))

#' GAM model with restricted ubic splines function fitted on the log2 and IQR transformed prenatal concentration of Cu
#'
#' @param sdq A string defining which SDQ score will be fitted
#' @param expo_list A character vector defining which exposures concentration will be fitted
#' @param input_data A dataframe with one imputed dataset
#' @param conf_list A character vector specifying names of model adjustment factors (confounders)
#' @param nknots A scalar defining number of knots for cubic splines
#'
#' @return Return a list of GAM models for each exposure of interest
#' @export
#'
#' @importFrom stats as.formula
#  @importFrom splines bs
#' @importFrom MASS negative.binomial
#' @importFrom gam gam
#' @importFrom rms rcs

GAMrcs <- function(input_data, conf_list, sdq, expo_list, nknots) {

  res <- list()

  for (expo in expo_list) {

  # Make a formula for the GAM equation (with cubic splines fitted on Cu)
  # formula_GAM <- stats::as.formula(paste(sdq, " ~ splines::bs(", expo, ", knots = 3) +", paste(conf_list, collapse = " + ")))
  formula_GAM_rcs <- stats::as.formula(paste(sdq, " ~ rms::rcs(", expo, ", knots = ", nknots, ") +", paste(conf_list, collapse = " + ")))

    # Run GAM
  # GAM <- gam::gam(formula = formula_GAM,
                  # family = negative.binomial(1),
                  # data = input_data)

  # Run GAM restricted cubic splines
  GAM_rcs <- gam::gam(formula = formula_GAM_rcs,
                  family = negative.binomial(1),
                  data = input_data)

  res[[length(res) + 1]] <- GAM_rcs

  }

  names(res) <- expo_list
  return(res)
}




