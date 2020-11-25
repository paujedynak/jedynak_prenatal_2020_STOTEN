# Manually set lists of compounds, models variables, etc. for Jedynak et al., 2020 paper
# Paulina Jedynak
# 23/09/20

# List of covariates names
covariates_list <- function() {
  covariates_list <- c("smoking", "mother_bmi", "mother_age", "mother_edu", "parity",
                       "breastf_cat", "mother_work", "child_age", "weight_gain",
                       "weight_gain_cat", "fish_intake", "conception_trim", "crowding",
                       "year_birth", "sex")

  return(covariates_list)
}


# List of variables that will be compared between cohorts
variables_to_compare_by_cohort <- function() {
  variables_to_compare_by_cohort <- c("cohort", "conception_trim", "smoking", "parity",
                                      "mother_edu", "mother_work", "mother_bmi", "weight_gain_cat",
                                      "sex", "child_age", "mother_age", "h_sdq_external",
                                      "h_sdq_internal")
  return(variables_to_compare_by_cohort)
}

# List of adjustment factors used in statistical analyses
conf_list <- function() {
  conf_list <- c("cohort", "smoking", "mother_bmi", "mother_age", "mother_edu", "parity",
                 "conception_trim", "mother_work", "child_age", "sex")
  return(conf_list)
}


# List of variables to be adjusted for fat
fat_adjust <- function() {
  fat_adjust <- c("dde", "ddt", "hcb", "pbde47", "pbde153", "pcb118", "pcb138", "pcb153",
                  "pcb170", "pcb180")

  return(fat_adjust)
}

# List of variables to be adjusted for creatinine
creatinine_adjust <- function() {
  creatinine_adjust <- c("mep", "mibp", "mnbp", "mbzp", "mehp", "mehhp", "meohp", "mecpp",
                         "ohminp", "oxominp", "mepa", 'etpa', "prpa", "bpa", "bupa", "oxbe",
                         "trcs", "dmp", "dmtp", "dmdtp", "dep", "detp", "cotinine")

  return(creatinine_adjust)
}

# List of variables to be log2 transformed
log2_transform <- function() {
  log2_transform <- c("pfoa", "pfna", "pfunda", "pfhxs", "pfos", "k", "mg", "na", "as", "cd",
                      "co", "cs", "cu", "hg", "mn", "mo", "pb", "se", "zn")
  return(log2_transform)
}


# List of variables to be ln transformed
ln_transform <- function() {
  ln_transform <- "hs_cotinine_adj"
  return(ln_transform)
}


# Load of original names of exposures along with short names of compounds and family name (ending with _m)
variable_list_m <- function() {
  variable_list_m <- readr::read_csv(file = here::here("data/variable_lists/variable_list_m.csv"))
  return(variable_list_m)
}

# List of exposure variables that will be used for regressions and other analysis, adjusted, log and IQR transformed (ending with _adj_log_iqr)
variable_list_adj_log_iqr <- function() {
  variable_list_adj_log_iqr <- readr::read_csv(here::here("data/variable_lists/variable_list_adj_log_iqr.csv"))
  return(variable_list_adj_log_iqr)
}

# List of exposure variables from OC family only
variable_list_oc <- function() {
  variable_list_oc <- dplyr::filter(variable_list_adj_log_iqr(), Variable_name %in%
                                      c("hs_dde_adj_log2_iqr", "hs_ddt_adj_log2_iqr",
                                        "hs_hcb_adj_log2_iqr", "hs_pcb118_adj_log2_iqr",
                                        "hs_pcb138_adj_log2_iqr", "hs_pcb153_adj_log2_iqr",
                                        "hs_pcb170_adj_log2_iqr", "hs_pcb180_adj_log2_iqr"))
  return(variable_list_oc)
}

# Variables used in multiple imputation
# The lists contain all the variables that will be treated specially during imputation (variables
#that will not be predictors, highly correlated variables, low predictability variables, skipped from
#imputation but not from prediction, creatinine adjusted, total fat percentage adjusted, not adjusted
#but log2 transformed)

# List variables that are highly correlated with each other
high_corr <- function() {
  high_corr <- c("pcb118", "pcb138", "pcb153", "pcb170", "pfos", "pfoa", "pfunda", "mnbp", "mibp", "mehhp", "meohp", "mecpp", "etpa", "mepa", "dmtp")
  return(high_corr)
}

high_corr_short <- function() {
  high_corr_short <- c("pfos", "pfoa", "pfunda")
  return(high_corr_short)
}


high_corr_dehp <- function() {
  high_corr_dehp <- c("dehp_sum_m", "dehp_sum_adj_log2")
  return(high_corr_dehp)
}

# List variables with low potential for prediction in imputation
low_pred <- function() {
  low_pred <- c("ddt", "pbde47", "pbde153", "as", "k", "mg", "na", "cd", "co", "cs", "cu", "mn", "mo", "pb", "se", "zn", "ohminp", "oxominp")
  return(low_pred)
}


low_pred_short <- function() {
  low_pred_short <- c("ddt", "pbde47", "pbde153", "ohminp", "oxominp")
  return(low_pred_short)
}

# List variables that will not be considered as predictors in imputation
no_pred <- function() {
  no_pred <- c("HelixID", "year_birth")
  return(no_pred)
}

# List variables to be skipped from imputation but not from prediction
skip <- function() {
  skip <- c("HelixID", "h_sdq_external", "h_sdq_internal", "cohort", "child_age", "sex", "year_birth")
  return(skip)
}

# Define which variables will be forced in the LASSO model (their coefficient won't be shrunk)
penalty_factor <- function() {
penalty_factor = rep(c(1, 0), times = c(47, 19))# 47: number of variables to be shrunk (exposures);
# 19: number of covariate variables (with contrasts added for categorical variables) that will not be
# shrunk
return(penalty_factor)
}
