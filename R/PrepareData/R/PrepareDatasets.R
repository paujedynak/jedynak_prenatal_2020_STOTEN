# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("h_cohort", "breastf_cat", "child_age", "cohort", "e3_asmokyn_p_None", "e3_sex_None",
                         "e3_yearbir_None", "h_age_None", "h_agemontsdq", "h_crowding_None", "h_edumc_None",
                         "h_fish_preg_Ter", "h_mbmic_None", "h_mwork_None", "h_parity_None", "h_sdq_external",
                         "h_sdq_internal", "h_trimcon_None", "hs_bf_None", "hs_creatinine_m", "hs_creatinine_mdesc",
                         "hs_creatinine_mg", "hs_ldlchol_m", "hs_ldlchol_mdesc", "hs_mecpp_m", "hs_mehhp_m",
                         "hs_mehp_m", "hs_meohp_m", "hs_phospho_m", "hs_phospho_mdesc", "hs_sumDEHPm_m", "hs_sumPCBs5_m",
                         "hs_totfat_mperc", "hs_wgtgain_None", "sex", "smoking", "weight_gain_cat", "mother_bmi",
                         "dehp_sum_m", "hs_k_m", "hs_mg_m", "hs_na_m", "hs_se_m", "hs_zn_m", "hs_k_mdesc", "hs_mg_mdesc",
                         "hs_na_mdesc", "hs_se_mdesc", "hs_zn_mdesc"))

#' Prepares all raw datasets for further processing
#' @param exposure_data A .csv file containing maternal exposures (data.frame, rows = subjects, cols = exposure variables)
#' @param covariate_data An .Rdata file containing covariates (data.frame, rows = subjects, cols = covariate variables)
#' @param sex_age_data A .csv file containing child age and sex covariates (data.frame, rows = subjects, cols = variables containing sex and age)
#' @param covariates_list A character vector containing covariates names
#' @param fat_adjust A character vector of names of variables to be adjusted with fat
#' @param creatinine_adjust A character vector of names of variables to be adjusted with creatinine
#' @param log2_transform A character vector of names of variables to be log2 transformed
#' @param ln_transform A character vector of names of variables to be ln transformed
#'
#' @return A dataframe with all outcomes, exposures and covariates of interest
#'
#' @export
#'
#' @import dplyr
#' @importFrom tibble add_column
#' @importFrom forcats fct_recode
#' @importFrom utils read.csv
#' @importFrom stats relevel
#' @importFrom tidyselect all_of
#'
PrepareDatasets <- function(exposure_data,
                            covariate_data,
                            covariates_list,
                            sex_age_data,
                            fat_adjust,
                            creatinine_adjust,
                            log2_transform,
                            ln_transform) {

  # Prepare the exposure-SDQ dataset for further analyses

  # Load and prepare dataset containing SDQ scores and biomarkers concentrations
  # For description of biomarkers refer to HELIX_Codebook Biomarkers_Mother_EXTERNAL_20170720

  SDQ_biomarkers <- exposure_data %>%

    # Select variables: cohort, HelixID, SDQ scores, total fat percentage, creatinine (g/L),
    # _m (pg/g serum in mother, numerical, Values below the limit of detection imputed),
    # _mraw and (raw in mother, numerical, only Tl and DEDTP, Raw data sent by the laboratory - values below LOD not imputed)
    # _mdesc (Sample description for each mother, categorical, Concentration above the limit of detection)

    dplyr::select(cohort, HelixID, h_sdq_external, h_sdq_internal,
                  hs_totfat_mperc, hs_creatinine_mg, ends_with("_m"),
                  ends_with("_mdesc"), ends_with("_mraw")) %>%

    # Recode SDQ scores to int
    dplyr::mutate(h_sdq_external = as.integer(h_sdq_external),
                  h_sdq_internal = as.integer(h_sdq_internal)) %>%

    # Change SAB cohort name to INMA
    dplyr::mutate(cohort = forcats::fct_recode(cohort, "INMA" = "SAB"),
                  cohort = factor(cohort,
                                  levels = c("BIB", "EDEN", "INMA", "KANC", "RHEA"),
                                  labels = c("BiB", "EDEN", "INMA", "KANC", "RHEA"))) %>%

    # Remove sumPCBs, sumDEHP, hs_creatinine_m, hs_phospho_m,
    # hs_totchol_m, hs_triglyc_m, hs_hdlchol_m and hs_ldlchol_m,
    # hs_phospho_mdesc, hs_totchol_mdesc, hs_triglyc_mdesc,
    # hs_hdlchol_mdesc and hs_ldlchol_mdesc variables,
    # hs_creatinine_mdesc, as these variables are not needed
    dplyr::select(-c(hs_sumPCBs5_m,
                     hs_sumDEHPm_m,
                     hs_creatinine_m,
                     hs_phospho_m:hs_ldlchol_m,
                     hs_phospho_mdesc:hs_ldlchol_mdesc,
                     hs_creatinine_mdesc)) %>%

    # Change names of creatinine and fat
    dplyr::rename(fat = hs_totfat_mperc, creatinine = hs_creatinine_mg) %>%

    # Create molar sum DEHP and add it after the DEHP metabolites -
    # we won't use the DEHP sum in the analyses but it will be used for imputation
    tibble::add_column(dehp_sum_m = NA, .after = "hs_mecpp_m") %>%
    dplyr::mutate(dehp_sum_m = hs_mehp_m/278 + hs_mehhp_m/294 + hs_meohp_m/292 + hs_mecpp_m/308)

  # Returns a dataset with 753 subjects (with missing values for SDQ scores)
  # and HelixID, cohort, 2 SDQ scores (internal and external), 2 adjustment
  # variables (fat and creatinine), 54 maternal exposures in numeric (_m or _mraw)
  # form and 54 maternal exposures in descriptive (_mdesc) form and
  # dehp sum variable (115 variables in total); N = 753


  # Prepare covariates dataset for further analyses

  # Drop imputed observations (leave only imp = 0 which is original dataset)
  covar <- dplyr::filter(covariate_data, .imp == 0) %>% # n = 1301
    dplyr::mutate(cohort = factor(h_cohort,
                                  levels = c("BIB", "EDEN", "INMA", "KANC", "MOBA", "RHEA"),
                                  labels = c("BiB", "EDEN", "INMA", "KANC", "MoBa", "RHEA")))

  cat("Women enrolled in HELIX study (total n): ")
  cat(dim(covar)[1])

  cat("\nWomen enrolled in HELIX study (n per cohort):")
  print(table(covar$cohort))


  cov <- covar %>%
    # Merge covariate_data and sex_age_data files, by 'HelixID'
    merge(sex_age_data, by = "HelixID") %>% # n = 746

    # Change variable names and correct errors, code factor variables
    dplyr::mutate(

      # Categorical variables
      smoking = factor(e3_asmokyn_p_None, labels = c("No", "Yes")),
      mother_edu = factor(h_edumc_None, labels = c("Primary", "Secondary", ">=University")),
      mother_work = factor(h_mwork_None, labels = c("Unemployed", "Employed")),
      mother_bmi = factor(h_mbmic_None, labels = c("Underweight", "Normal_weight", "Overweight", "Obesity")),
      mother_bmi = stats::relevel(mother_bmi, ref = "Normal_weight"),
      parity = factor(h_parity_None, labels = c("No_child", "1_child", "2_or_more_children")),
      conception_trim = factor(h_trimcon_None, labels = c("Jan-March", "April-June", "July-Sept", "Oct-Dec")),
      fish_intake = factor(h_fish_preg_Ter, labels = c("T1", "T2", "T3")), #T1: <1.9 / T2: 1.9-4.1 / T3: >4.1
      sex = factor(e3_sex_None, labels = c("Female", "Male")),
      year_birth = factor(e3_yearbir_None),

      # Numerical variables
      child_age = h_agemontsdq / 12,
      mother_age = h_age_None,
      weight_gain = hs_wgtgain_None,
      crowding = h_crowding_None,

      # Correct mistakes spotted in breastfeeding variable
      breastf_cat = hs_bf_None,
      breastf_cat = dplyr::case_when(
        HelixID == "SAB60" |
          HelixID == "EDE7496" |
          HelixID == "EDE4269" |
          HelixID == "EDE7150" |
          HelixID == "SAB646" ~ "1",

        HelixID == "EDE4070" ~ "0",

        TRUE ~ as.character(hs_bf_None)),
      breastf_cat = factor(breastf_cat, labels = c("No", "Yes")),

      # Create categorical weight gain dependent on pre-pregnancy BMI variable
      weight_gain_cat = dplyr::case_when(
        mother_bmi == "Underweight" & weight_gain < 12.5 |
          mother_bmi == "Normal_weight" & weight_gain < 11.5 |
          mother_bmi == "Overweight" & weight_gain < 7 |
          mother_bmi == "Obesity" & weight_gain < 5 ~ "Insufficient",

        mother_bmi == "Underweight" & weight_gain > 18 |
          mother_bmi == "Normal_weight" & weight_gain > 16 |
          mother_bmi == "Overweight" & weight_gain > 11.5 |
          mother_bmi == "Obesity" & weight_gain > 9 ~ "Excessive",

        TRUE & !is.na(weight_gain) & !is.na(mother_bmi) ~ "Adequate"),
      weight_gain_cat = factor(weight_gain_cat, levels = c("Insufficient", "Adequate", "Excessive")),
      weight_gain_cat = stats::relevel(weight_gain_cat, ref = "Adequate")
    )

  #=================

  # Select previously chosen covariate variables

  # Select covariates (excl. cohort) and HelixID for further analysis (16 variables, N = 746)
  covariates <- cov %>%
    dplyr::select(HelixID, tidyselect::all_of(covariates_list))

  # Merge datasets to obtain final dataset that will be used in further analyses
  # Merge biomarkers-SDQ dataset and covariates by HelixID
  sdq_expo_cov <- SDQ_biomarkers %>%
    dplyr::left_join(covariates, by = "HelixID") %>%

    # Remove incomplete cases from sdq internal, sdq external, sex and age
    dplyr::filter(stats::complete.cases(h_sdq_external) &
                    stats::complete.cases(h_sdq_internal)) # n = 721

  cat("\nSDQ scores available (total n): ")
  cat(dim(sdq_expo_cov)[1])

  cat("\nSDQ scores available (n per cohort):")
  print(table(sdq_expo_cov$cohort))

  sdq_expo_cov <- sdq_expo_cov %>%
    dplyr::filter(stats::complete.cases(sex) &
                    stats::complete.cases(child_age)) %>%

    # Returns a dataset with 708 subjects (final n. Complete cases for SDQs, age and sex)
    # and HelixID, cohort, 2 SDQ scores, 2 adjustment variables (fat and creatinine),
    # 52 maternal exposures in numeric (_m or _mraw) form, sum_DEHP and 15 covariates
    # (74 variables in total)

    # Prepare a final dataset that will not contain Tl or DEDTP and descriptive ("_mdesc" variables)
    # Remove hs_tl_mraw and hs_dedtp_mraw as missing for most of the variables (we won't take Tl and
    # DEDTP into account in further analyses). Remove _mdesc variables as not needed either
    dplyr::select(-contains("tl"), -contains("dedtp"), -ends_with('_mdesc')) %>%
    as.data.frame()

  cat("\nStudy population (total n): ")
  cat(dim(sdq_expo_cov)[1])

  cat("\nStudy population (n per cohort):")
  print(table(sdq_expo_cov$cohort))

  # For most of the analyses, exposure variables need to be adjusted for fat/ creatinine and log2/ln transformed (where appropriate).
  # Adjust exposure concentrations for fat or creatinine, where appropriate

  # Make a list of variables to be adjusted for fat
  fat_list <- ChangeVarsNames(x = fat_adjust, hs_m = TRUE)

  # Make a list of variables to be adjusted for creatinine
  creatinine_list <- ChangeVarsNames(x = creatinine_adjust, hs_m = TRUE) %>%
    append("dehp_sum_m")

  sdq_expo_adj_cov <- .AdjustFatCreat(data = sdq_expo_cov,
                                      list_fat_adj = fat_list,
                                      list_creat_adj = creatinine_list)

  #=====================

  # Log2 or ln exposure concentrations, where appropriate

  # Make a list of variables to be log2 transformed
  log2_list <- ChangeVarsNames(x = log2_transform, hs_m = TRUE)

  # There is no need to manipulate this list as there is only one variable
  ln_list <- ln_transform

  sdq_expo_adj_log_cov <- .LogTransform(data = sdq_expo_adj_cov,
                                                    list_log2 = log2_list,
                                                    list_ln = ln_list)

  # From the SDQ_biomarkers dataset remove the non-toxic elements (K, Mg, Na, Se, Zn) and DEHP sum.
  # They will be kept in the transformed dataset as they will be used for imputation

  SDQ_biomarkers <- SDQ_biomarkers %>%
    dplyr::select(cohort, HelixID, ends_with('_m'), ends_with('_mraw'), ends_with('_mdesc'),
                  -c(fat, creatinine, hs_k_m, hs_mg_m, hs_na_m, hs_se_m, hs_zn_m, hs_k_mdesc, hs_mg_mdesc, hs_na_mdesc, hs_se_mdesc, hs_zn_mdesc, dehp_sum_m)) %>%

    # Select complete cases in regard to SDQ scores, child sex an age
    dplyr::filter(HelixID %in% sdq_expo_cov$HelixID)


  # Create a dataset with exposure descriptive variables only
  exposures_description <- SDQ_biomarkers %>%
    dplyr::select(ends_with("_mdesc"))

  # Create a dataset with exposure concentrations variables only
  exposures_concentrations <- SDQ_biomarkers %>%
    dplyr::select(cohort, ends_with("_m"))

  # From the sdq_expo_cov dataset remove the non-toxic elements (K, Mg, Na, Se, Zn) and DEHP sum.
  # They will be kept in the transformed dataset as they will be used for imputation
  sdq_expo_cov <- sdq_expo_cov %>%
    dplyr::select(-c(dehp_sum_m, hs_k_m, hs_mg_m, hs_na_m, hs_se_m, hs_zn_m))

  # Save an object with all datasets of interest
  recoded_datasets <- list("exposures_description" = exposures_description,
                           "exposures_concentrations" = exposures_concentrations,
                           "sdq_expo_cov" = sdq_expo_cov,
                           "sdq_expo_adj_log_cov" = sdq_expo_adj_log_cov)

  # Save dataset into global environment
  return(recoded_datasets)
}
