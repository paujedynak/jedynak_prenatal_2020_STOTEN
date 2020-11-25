# Tables for Jedynak et al., 2020 paper
# Paulina Jedynak
# 23/09/20

# Table 1
Tab_1 <- function(res,
                    res_coh,
                    path,
                    file_name) {

  Table_1 <- cbind(res, res_coh)
  colnames(Table_1) <- c("Covariate", "Overall", "EDEN", "BiB", "INMA", "KANC", "RHEA", "p_value")

  Table_1 <- Table_1 %>%
    dplyr::select(Covariate, Overall, BiB, EDEN, INMA:p_value)

  utils::write.csv(Table_1, here::here(path, paste0(file_name, ".csv")))

  return(Table_1)
}

# Table 2
Tab_2 <- function(exwas_ext,
                    exwas_int,
                    MLR_ext,
                    MLR_int,
                    complete_case_ext,
                    complete_case_int,
                    expo_sdq_ext,
                    expo_sdq_int,
                    exwas_ext_no_BiB,
                    n_no_BiB,
                    path,
                    file_name) {

  main_ExWAS <- Helpers::MergeDatasets(input_data_ext = exwas_ext,
                                      input_data_int = exwas_int,
                                      expo_sdq_ext = expo_sdq_ext,
                                      expo_sdq_int = expo_sdq_int)

  MLR <- Helpers::MergeDatasets(input_data_ext = MLR_ext,
                               input_data_int = MLR_int,
                               expo_sdq_ext = expo_sdq_ext,
                               expo_sdq_int = expo_sdq_int)

  compl_case <- Helpers::MergeDatasets(input_data_ext = complete_case_ext,
                                      input_data_int = complete_case_int,
                                      expo_sdq_ext = expo_sdq_ext,
                                      expo_sdq_int = expo_sdq_int)

  no_BiB = dplyr::select(exwas_ext_no_BiB, -c(Family, Exposure, p_value_FDR))

  Table_2 <- dplyr::left_join(main_ExWAS, MLR, by = "Variable_name") %>%
    dplyr::left_join(compl_case, by = "Variable_name") %>%
    dplyr::full_join(no_BiB, by = "Variable_name") %>%
    dplyr::mutate(n_no_BiB = n_no_BiB) %>%
    dplyr::rename(estCI_ExWAS = estCI.x,
                  p_value_ExWAS = p_value.x,
                  p_value_FDR_ExWAS = p_value_FDR,
                  estCI_MLR = estCI.y,
                  p_value_MLR = p_value.y,
                  estCI_complete_c = estCI.x.x,
                  p_value_complete_c = p_value.x.x,
                  n_obs_complete_c = n_obs,
                  estCI_no_BiB = estCI.y.y,
                  p_value_no_BiB = p_value.y.y) %>%
    dplyr::mutate(SDQ = dplyr::case_when(Variable_name %in% expo_sdq_ext ~ "SDQ_external",
                                         Variable_name %in% expo_sdq_int ~ "SDQ_internal")) %>%
    dplyr::select(SDQ, Exposure, Family, estCI_ExWAS:n_no_BiB) %>%
    dplyr::arrange(SDQ, Exposure)

  utils::write.csv(Table_2, here::here(path, paste0(file_name, ".csv")))

  return(Table_2)
}

# Appendix Table 5
App_Tab_5 <- function(exposure_char,
                      variable_list_m,
                      path,
                      file_name) {

  App_Tab_5 <- dplyr::left_join(exposure_char$miss_LOD, exposure_char$median_IQR_overall, by = "Variable_short") %>%
    dplyr::left_join(exposure_char$median_IQR_per_cohort, by = "Variable_short") %>%
    dplyr::left_join(exposure_char$`p_value_Kruskall-Wallis`, by = "Variable_short") %>%
    dplyr::right_join(variable_list_m, ., by = "Variable_short") %>%
    dplyr::arrange(Family, Exposure) %>%
    dplyr::select(-c(Variable_name, Variable_short, Family))

  utils::write.csv(App_Tab_5, file = here::here(path, paste0(file_name, ".csv")))

  return(App_Tab_5)
}

# Appendix Tables 6 & 7
App_Tab_6_7 <- function(exwas,
                        heterog_estimates,
                        path,
                        file_name) {

  App_Tab_6_7 <- merge(dplyr::select(exwas, -Variable_name),
                       dplyr::select(heterog_estimates, -Family), by = "Exposure")

  utils::write.csv(App_Tab_6_7, file = here::here(path, paste0(file_name, ".csv")))

  return(App_Tab_6_7)
}


