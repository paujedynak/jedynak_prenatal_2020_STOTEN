# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("GVIF",
                         "Variable_name",
                         "term",
                         "conf_low",
                         "conf_high",
                         "p_value",
                         "estCI"))

#' Calculate Variance inflation factor (VIF) for a given model
#'
#' @param input_data A dataframe with a complete dataset
#' @param sdq A string indicating which SDQ score will be plugged into the model
#' @param var_list A character vector containing a list of exposures for which regression will be run
#' @param sign_var_list A character vector containing a list of exposures detected as significanly associated with SDQ in ExWAS
#' @param conf_list A character vector specifying confounders to be added in the regression model
#'
#' @return A character vector with exposures that are not collinear and can be plugged into a multpile regression model
#' @export
#' @import dplyr
#' @importFrom tidyselect all_of
#' @importFrom stats as.formula
#' @importFrom car vif
#' @importFrom MASS glm.nb
#' @importFrom Hmisc %nin%
#'
CalculateVIF <- function(input_data, sdq, var_list, sign_var_list, conf_list) {

repeat {

    # Select variables list to be plugged in the regression formula
    formula_VIF <- stats::as.formula(paste(sdq, "~", paste(var_list, collapse = " + "), "+", paste(conf_list, collapse = " + ")))

    # Fit a negative binomial model
    model <- MASS::glm.nb(formula = formula_VIF, data = input_data)

    # Create a model matrix on which VIFs will be calculated
    vif <- car::vif(model) %>% # GVIFs over 4 indicate collinear variables
      as.data.frame() %>%

      # Filter out the confounder variables
      dplyr::filter(rownames(.) %nin% c(sign_var_list, conf_list))

    # List the exposure to be removed because of collinearity
    to_remove <- dplyr::filter(vif, GVIF >= 4 & GVIF == max(GVIF))
    remove_var <- rownames(to_remove)

    if (length(remove_var) != 0) {

      print(paste0("Variable with high VIF: ", remove_var))
      print(paste0("VIF: ", to_remove$GVIF))
      var_list <- var_list[var_list != remove_var]

    } else {

     break

    }

  }

  return(var_list)

}
