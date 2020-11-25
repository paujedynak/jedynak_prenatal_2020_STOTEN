# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("Concentration", "weight_gain_cat"))

#' Calculate regression estimate (ExWAS) for each category of data stratified on weight gain
#'
#' @param data A data frame of complete dataset
#' @param strat_var_list A character vector containing the names of variables to stratify on
#' @param variable_list A data.frame containing a list of exposures for which regressions will be run, first column must contain exposure names
#'
#' @return A dataframe with regression estimates with p values for each exposure of interest
#' @export

StratifyGainWeight <- function(data, strat_var_list, variable_list) {

  # Make a long dataset including all the variables of interest: weight gain based on the BMI, OC family members and sdq external
res <- data %>%
  tidyr::gather(key = Variable_name, value = Concentration, strat_var_list$Variable_name) %>%
  tidyr::gather(key = weight_gain_cat, value = group, weight_gain_cat) %>%

  # Run a stratified analysis according to the weight gain in reference to pre-pregnancy BMI (group) for each member of the OC family
  dplyr::group_by(group, Variable_name) %>%
  dplyr::group_modify(~ broom::tidy(MASS::glm.nb(h_sdq_external ~
                                                   Concentration + mother_bmi + cohort + smoking + mother_age + mother_edu + parity + conception_trim + mother_work + child_age + sex, data = .),
                                    conf.int = TRUE,
                                    exponentiate = TRUE)) %>%
  dplyr::filter(term == "Concentration") %>%
  merge(dplyr::select(variable_list, Variable_name, Exposure), by = "Variable_name", .) %>%

  # Select variables to display
  dplyr::select(group, Variable_name, Exposure, estimate, conf.low, conf.high, p.value) %>%
  dplyr::rename(conf_low = conf.low, conf_high = conf.high, Estimate = estimate, p_value = p.value) %>%
  dplyr::arrange(Variable_name)

return(res)
}
