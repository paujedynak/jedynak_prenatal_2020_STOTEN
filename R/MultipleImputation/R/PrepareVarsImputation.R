#' Changes names of the exposures that will be used in the imputation
#'
#' @param high_corr A character vector of highly correlated variables
#' @param high_corr_short A character vector with a subset of highly correlated variables
#' @param high_corr_dehp A character vector with a subset of highly correlated variables from DEHP sum group
#' @param low_pred A character vector of variables with low predictability
#' @param low_pred_short A character vector with a subset of variables with low predictability
#' @param creatinine_adjust A character vector of variables adjusted for creatinine
#' @param fat_adjust A character vector of variables adjusted for fat
#' @param log2_transform A character vector of log transformed variables
#'
#' @return A list whose each element contains a list of newly named variables
#' @export
#' @importFrom Hmisc %nin%
#' @importFrom PrepareData ChangeVarsNames
#'
#'
PrepareVarsImputation <- function(high_corr,
                                  high_corr_short,
                                  high_corr_dehp,
                                  low_pred,
                                  low_pred_short,
                                  creatinine_adjust,
                                  fat_adjust,
                                  log2_transform) {

  # List highly correlated variables
  high_corr_list <- sapply(high_corr, FUN = PrepareData::ChangeVarsNames, hs_m = TRUE, USE.NAMES = FALSE) %>%
    append(sapply(high_corr[high_corr %nin% high_corr_short], PrepareData::ChangeVarsNames, adj_log2 = TRUE, USE.NAMES = FALSE)) %>%
    append(sapply(high_corr_short, PrepareData::ChangeVarsNames, m_log2 = TRUE, USE.NAMES = FALSE)) %>%
    append(high_corr_dehp)

  # List variables with low potential for prediction in imputation
  low_pred_list <- sapply(low_pred, PrepareData::ChangeVarsNames, hs_m = TRUE, USE.NAMES = FALSE) %>%
    append(sapply(low_pred[low_pred %nin% low_pred_short], PrepareData::ChangeVarsNames, m_log2 = TRUE, USE.NAMES = FALSE)) %>%
    append(sapply(low_pred_short, function(x) paste0("hs_", x, "_adj_log2"), USE.NAMES = FALSE))

  # List creatinine adjusted variables
  creat_list <- creatinine_adjust[creatinine_adjust != "cotinine"]
  creat_list <- creat_list %>%
    lapply(function(x) c(paste0("hs_", x, "_adj_log2"), paste0("~ I(log2(hs_", x, "_m / creatinine))"), paste0("hs_", x, "_m"), paste0("hs_", x, "_adj"))) %>%

    # add variables that have different name pattern (log2 vs ln)
    append(list(c("dehp_sum_adj_log2", "~ I(log2(dehp_sum_m / creatinine))", "dehp_sum_m", "dehp_sum_adj")), after = 10) %>%
    append(list(c("hs_cotinine_adj_ln", "~ I(log(hs_cotinine_m / creatinine))", "hs_cotinine_m", "hs_cotinine_adj")), after = length(.))

  # List total fat adjusted variables
  fat_list <- fat_adjust %>%
    lapply(function(x) c(paste0("hs_", x, "_adj_log2"), paste0("~ I(log2((hs_", x, "_m / 2) / (fat * 10)))"), paste0("hs_", x, "_m"), paste0("hs_", x, "_adj")))

  # List variables that are not fat/creat adjusted but log2 transformed (PFASs and metals)
  log2_list <- log2_transform %>%
    lapply(function(x) c(paste0("hs_", x, "_m_log2"), paste0("~ I(log2(hs_", x, "_m))"), paste0("hs_", x, "_m")))

  vars_imput <- list("high_corr_list" = high_corr_list,
                     "low_pred_list" = low_pred_list,
                     "creat_list" = creat_list,
                     "fat_list" = fat_list,
                     "log2_list" = log2_list)
  return(vars_imput)
  }

