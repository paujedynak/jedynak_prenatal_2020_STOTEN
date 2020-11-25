# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("SDQ",
                         "exposure",
                         "term",
                         "p.value"))

#' helper function to calculate interaction between sex and exposures identified in the main ExWAS
#'
#' @param input_data A mids object (imputed dataset)
#' @param sdq A string characterizing which SDQ score is to be used "h_sdq_external" or "h_sdq_internal"
#' @param conf_list A character vector listing the adjustment factors
#' @param strat_var_list A character vector listing the names of variables for which interaction with sex will be checked
#'
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom mice complete
#' @importFrom stats as.formula
#' @importFrom broom tidy
#'
#' @return A dataframe with p value for interaction between sex and exposure of interest in an ExWAS adjusted for confounding factors
#' @export

InteractionSex <- function(input_data, sdq, conf_list, strat_var_list) {

  # Create a formula for regression (as sex will be in the interaction term, remove it from the rest of the formula)
  formula <- stats::as.formula(paste("score ~ concentration * sex", "+", paste(conf_list[conf_list != "sex"], collapse = " + ")))

  sex_interaction <- mice::complete(input_data) %>%
    tidyr::gather(key = "SDQ", value = "score", sdq) %>%
    tidyr::gather(key = "exposure", value = "concentration", strat_var_list) %>%
    dplyr::group_by(SDQ, exposure) %>%
    dplyr::group_modify(~ broom::tidy(MASS::glm.nb(formula = formula, data = .), conf.int = TRUE, exponentiate = TRUE)) %>%
    dplyr::filter(term == "concentration:sexMale") %>%
    dplyr::select(SDQ, exposure, p.value) %>%
    dplyr::arrange(SDQ, p.value) %>%
    dplyr::ungroup()

  return(sex_interaction)
}
