# Jedynak et al., 2020 paper analysis drake plan
# Paulina Jedynak
# 23/09/20

plan <- drake::drake_plan(

  # Data import ----

  # Load dataset containing exposures and SDQ scores data
  exposure_data = readr::read_csv(file = here::here("data/raw_data/SDQ_biomarkers_july_2018.csv")),

  # Load dataset containing covariates data
  covariates_data = readRDS(here::here("data/raw_data/20180212_v2_3_ImputedDataset_20x10_AllDataFrame.RDS")),

  # Load dataset containing sex and age of the children
  sex_age_data = utils::read.csv(file = here::here("data/raw_data/INMA_EDEN_KANK_RHEA_BIB_SDQ_2018_07.csv"), stringsAsFactors = TRUE),

  # Data management ----

  # Before further analyses, original datasets need to be pre-processed: variables of interest need to
  # be selected, some variables need to be renamed and recoded, missing values removed or imputed,
  # datasets merged etc.

  # We divided urinary biomarker concentrations by creatinine concentration. Haemal lipophilic biomarker
  # concentrations were standardized and expressed in ng/g of total lipids in serum or plasma.
  # Concentrations were then ln-transformed (cotinine) or log2-transformed (all other biomarkers) to
  # approach normality and standardized for the interquartile range (IQR) by dividing exposure
  # concentration observed for each individual for a given exposure by the IQR calculated for this exposure.

  # Merge all original datasets (exposures, SDQ scores, covariates) in one dataset.
  # This operation will limit the number of subjects to those with complete SDQ scores, child age and sex. It will rename covariates of interest and create weight gain and DEHP sum variables.
  sdq_expo_cov_datasets = PrepareData::PrepareDatasets(exposure_data = exposure_data,
                                                       covariate_data = covariates_data,
                                                       covariates_list = covariates_list(),
                                                       sex_age_data = sex_age_data,
                                                       fat_adjust = fat_adjust(),
                                                       creatinine_adjust = creatinine_adjust(),
                                                       log2_transform = log2_transform(),
                                                       ln_transform = ln_transform()),


  # 4 datasets will be created: exposures with LOD status (exposures_description), exposures with
  # concentration values (exposures_concentrations), combined dataset of covariates and SDQ and exposures
  # (sdq_expo_cov) and a combined dataset with SDQ, exposures fat/creat adjusted and log transformed and
  # covariates (sdq_expo_adj_log_cov)

  # Exposures with LOD status
  exposures_description = sdq_expo_cov_datasets$exposures_description,

  # Exposures with concentration values (exposures_concentrations)
  exposures_concentrations = sdq_expo_cov_datasets$exposures_concentrations,

  # Combined dataset of covariates and SDQ and exposures
  sdq_expo_cov = sdq_expo_cov_datasets$sdq_expo_cov,

  # Combined dataset with SDQ, exposures fat/creat adjusted and log transformed and covariates
  sdq_expo_adj_log_cov = sdq_expo_cov_datasets$sdq_expo_adj_log_cov,


  # Data imputation ----

  # Missing data for exposure biomarker concentrations and adjustment factors were multiply
  # imputed (100 imputed datasets) via a chained equations algorithm

  # Before the imputation, lists of variables need to be created that will contain variables' names
  # that will go through passive imputation or will be removed from predictors (e.g. because of high
  # correlation etc.)

  vars_imput = MultipleImputation::PrepareVarsImputation(high_corr = high_corr(),
                                                         high_corr_short = high_corr_short(),
                                                         high_corr_dehp = high_corr_dehp(),
                                                         low_pred = low_pred(),
                                                         low_pred_short = low_pred_short(),
                                                         creatinine_adjust = creatinine_adjust(),
                                                         fat_adjust = fat_adjust(),
                                                         log2_transform = log2_transform()),
  # Impute 100 datasets
  imputed_datasets = MultipleImputation::MultipleImputation(dataset = sdq_expo_adj_log_cov,
                                                            no_pred_list = no_pred(),
                                                            high_corr_list = vars_imput$high_corr_list,
                                                            low_pred_list = vars_imput$low_pred_list,
                                                            skip_list = skip(),
                                                            creat_list = vars_imput$creat_list,
                                                            fat_list = vars_imput$fat_list,
                                                            log2_list = vars_imput$log2_list,
                                                            m = 100,
                                                            seed = 1083024,
                                                            path = "results/main_analyses",
                                                            file_name = "imputed_datasets.RDS"),

  # Select a dataset with IQR transformed biomarker exposure concentrations
  imp_100_adj_log_iqr = imputed_datasets$imputed_IQR_data,


  # Create different versions of imputed dataset for further analyses

  # Create a stacked dataset from the imputed dataset (stack all imputed datasets together with
  # the original  dataset with missing values)
  stacked_dataset = Helpers::Mids2Stacked(input_data = imp_100_adj_log_iqr,
                                          variable_list = variable_list_adj_log_iqr(),
                                          include = TRUE) %>%
    Helpers::ChangeCohortOrder(),

  # Create a stacked dataset from the imputed dataset (stack all imputed but do not include original
  # dataset with missing values)
  stacked_dataset_imp_only = Helpers::Mids2Stacked(input_data = imp_100_adj_log_iqr,
                                                   variable_list = variable_list_adj_log_iqr(),
                                                   include = FALSE),

  # Create one complete dataset
  complete_dataset = Helpers::Mids2SingleComplete(input_data = imp_100_adj_log_iqr,
                                                  variable_list = variable_list_adj_log_iqr()),


  # Descriptive analysis ----

  # Population characteristics for the mother-child pairs included in the study: overall and by cohort

  population_char = DescriptiveStats::PopulationCharacteristics(
    input_data = sdq_expo_cov,
    variable_list = variables_to_compare_by_cohort(),
    path = "results/main_analyses"),

  # Maternal exposure descriptive statistics: overall and by cohort

  exposure_char = DescriptiveStats::ExposureCharacteristics(input_data_descr = exposures_description,
                                                            input_data_conc = exposures_concentrations),


  # Main statistical analyses ----

  # LASSO

  # We used a least absolute shrinkage and selection operator (LASSO) algorithm with log link function.
  # LASSO considers all exposures simultaneously [@tibshirani1996] and performs variable selection
  # through estimates’ shrinkage (i.e. the lowest regression coefficients corresponding to the least
  # informative predictors are assigned a zero value). We determined the range of penalty parameter λ by
  # maximizing the prediction log-likelihood using 10-fold cross-validation. To prevent overfitting, we
  # defined the optimal λ as the one providing the sparsest model (as measured by the number of nonzero
  # regression coefficients) among those yielding a log-likelihood within 1 standard error of the maximum
  # log-likelihood [@krstajic2014]. To stabilise estimates, LASSO was fit on each of the 100 imputed
  # datasets and an exposure was retained only if it was selected in at least 50% of runs [@wood2008].


  # SDQ external
  lasso_ext = LassoENET::Lasso(input_data = stacked_dataset_imp_only,
                               sdq = "h_sdq_external",
                               variable_list = variable_list_adj_log_iqr(),
                               conf_list = conf_list(),
                               penalty_factor = penalty_factor(),
                               path = "results/main_analyses",
                               file_name = "lasso_ext"),

  # SDQ internal
  lasso_int = LassoENET::Lasso(input_data = stacked_dataset_imp_only,
                               sdq = "h_sdq_internal",
                               variable_list = variable_list_adj_log_iqr(),
                               conf_list = conf_list(),
                               penalty_factor = penalty_factor(),
                               path = "results/main_analyses",
                               file_name = "lasso_int"),

  # Display which variables exceeded the 50% selection threshold

  # SDQ external
  hits_lasso_ext = LassoENET::HitsNumber(regr_result = lasso_ext,
                                         variable_list = variable_list_adj_log_iqr()),

  # SDQ internal
  hits_lasso_int = LassoENET::HitsNumber(regr_result = lasso_int,
                                         variable_list = variable_list_adj_log_iqr()),


  # Exposome-wide association study (ExWAS) adjusted for confounding factors

  # To compare with previous single-pollutant studies, we also performed an exposome-wide
  # association study (ExWAS): we fit a negative binomial regression model on each of the 100
  # imputed datasets for each exposure and SDQ score, then aggregated the results using Rubin’s
  # rule for multiply imputed data (Patel et al. 2010)⁠. To control for multiple comparisons, we
  # applied a family wise error rate (FWER) correction to the p value threshold. The correction
  # uses a Bonferroni procedure extended to handle correlated tests: the actual number of
  # exposures being tested (M) is replaced by a smaller value called the effective number of
  # independent exposures (Me). Me is estimated by ∑_(i=1)^M[I(λ_i>1)(λ_i-1)], where I(x) is
  # an indicator function and λ_i are the eigenvalues of the matrix of correlations between M
  # exposures. The p value threshold to control FWER to α, using Me in a Bonferroni procedure, is then
  # α / Me (adapted from [@li2012]).

  # SDQ external
  exwas_ext = ExWAS::RunExwasImputed(input_data = imp_100_adj_log_iqr,
                                     variable_list = variable_list_adj_log_iqr(),
                                     sdq = "h_sdq_external",
                                     conf_list = conf_list()),

  # SDQ internal
  exwas_int = ExWAS::RunExwasImputed(input_data = imp_100_adj_log_iqr,
                                     variable_list = variable_list_adj_log_iqr(),
                                     sdq = "h_sdq_internal",
                                     conf_list = conf_list()),

  # List exposures significantly associated with SDQ score

  # SDQ external
  expo_sdq_ext = ExWAS::SelectSignificant(data = exwas_ext,
                                          pvalue = 0.05),

  variable_list_ExWAS_ext = dplyr::filter(variable_list_adj_log_iqr(), Variable_name %in% expo_sdq_ext) %>%
    dplyr::arrange(Exposure),


  # SDQ internal
  # Note: pvalue = 0.055 because we want to include the exposures that were on the verge of significance
  expo_sdq_int = ExWAS::SelectSignificant(data = exwas_int,
                                          pvalue = 0.055),

  variable_list_ExWAS_int = dplyr::filter(variable_list_adj_log_iqr(), Variable_name %in% expo_sdq_int) %>%
    dplyr::arrange(Exposure),

  # We evaluated the between-cohort heterogeneity of the adjusted association using the I².
  # Cohort-specific interactions were accounted for by computing the I² statistic [@higgins2002].
  # The I² statistic evaluates the between-cohort heterogeneity of the association between each
  # exposure and SDQ externalizing and internalizing scores, adjusted for confounders. The I² was
  # computed by adding an interaction term between each exposure and the cohort origin in an adjusted
  # negative binomial regression model. Low I² values suggest low heterogeneity across cohorts
  # [@guyatt_grade_2011]. We relied on the following threshold for I2 interpretation: I2 < 0.3
  # for low heterogeneity, 0.3 ≤ I2 < 0.6 for moderate heterogeneity, I2 ≥ 0.6 for substantial
  # to high heterogeneity [@deeks2019].

  # Note: Code for calculation of the I² in the following section (Heterogeneity package) was
  # adapted from Xavier Basagaña.

  # SDQ external
  heterog_estimates_ext = Heterogeneity::CaluclateHeterogeneityEstimates(data = stacked_dataset,
                                                                         outcome = "h_sdq_external",
                                                                         conf_list = conf_list(),
                                                                         variable_list = variable_list_adj_log_iqr(),
                                                                         strat_variable = "cohort"),
  # SDQ internal
  heterog_estimates_int = Heterogeneity::CaluclateHeterogeneityEstimates(data = stacked_dataset,
                                                                         outcome = "h_sdq_internal",
                                                                         conf_list = conf_list(),
                                                                         variable_list = variable_list_adj_log_iqr(),
                                                                         strat_variable = "cohort"),


  # Sensitivity analyses ----

  # To test the robustness of the associations between SDQ scores and exposures identified by the
  # LASSO (selected in at least 50% of runs) and ExWAS (those whose uncorrected p values < 0.05) we
  # performed further sensitivity analyses.

  # We evaluated the linearity of the associations using generalized additive model (GAM) with
  # restricted cubic splines function.

  # Fit a GAM model for SDQ external
  GAM_ext = GAMrcs::GAMrcs(input_data = complete_dataset,
                           conf_list = conf_list(),
                           sdq = "h_sdq_external",
                           expo_list = expo_sdq_ext,
                           nknots = 3),

  # Fit a GAM model for SDQ external
  GAM_int = GAMrcs::GAMrcs(input_data = complete_dataset,
                           conf_list = conf_list(),
                           sdq = "h_sdq_internal",
                           expo_list = expo_sdq_int,
                           nknots = 3),

  # We ran a regression simultaneously adjusted for all exposures associated with the SDQ scores
  # in the main ExWAS (p values < 0.2).

  # SDQ external
  var_list_ext = ExWAS::SelectSignificant(data = exwas_ext,
                                          pvalue = 0.2),

  # Main analysis: BPA, Cu, MnBP, PCB-138, DDT, PFUnDA plus Cd, Co, DDE, HCB, PBDE-47, PCB-118, PCB-153,
  # PRPA

  # SDQ internal
  var_list_int = ExWAS::SelectSignificant(data = exwas_int,
                                          pvalue = 0.2),

  # Main analysis: DETP, PFOS plus BUPA, Co, Mn, PFHXS, PFNA, PFUnDA


  # We need to check the Variance inflation factor (VIF) for the exposures that were selected as
  # significantly associated with the SDQ score in ExWAS (p value < 0.2) to see if they all can be
  # plugged into multiple model (to avoid collinearity).

  # SDQ external
  mlr_var_list_ext = MultipleRegression::CalculateVIF(input_data = complete_dataset,
                                                      sdq = "h_sdq_external",
                                                      var_list = var_list_ext,
                                                      sign_var_list = expo_sdq_ext,
                                                      conf_list = conf_list()),

  # [1] "Variable with high VIF: hs_pcb153_adj_log2_iqr"
  # [1] "VIF: 18.72"
  #
  # Because the PCBs were very highly correlated with each other, PCB-153 will be removed basing on
  # VIF (> 4 means collinearity) and correlation coefficient to avoid collinearity.

  # SDQ internal
  mlr_var_list_int_temp = MultipleRegression::CalculateVIF(input_data = complete_dataset,
                                                                         sdq = "h_sdq_internal",
                                                                         var_list = var_list_int,
                                                                         sign_var_list = expo_sdq_int,
                                                                         conf_list = conf_list()),
  # [1] "Variable with high VIF: hs_pfna_m_log2_iqr"
  # [1] "VIF: 4.19"
  #
  # For PFNA the VIF was only slightly higher than 4 and this compound was not strongly correlated with
  # other compounds, so it will not be removed from the analyses

  mlr_var_list_int = var_list_int,

  # Run multiple models

  # SDQ external
  MLR_ext = MultipleRegression::MultipleRegression(input_data = imp_100_adj_log_iqr,
                                                   mlr_var_list = mlr_var_list_ext,
                                                   variable_list = variable_list_adj_log_iqr(),
                                                   conf_list = conf_list(),
                                                   sdq = "h_sdq_external"),

  # SDQ internal
  MLR_int = MultipleRegression::MultipleRegression(input_data = imp_100_adj_log_iqr,
                                                   mlr_var_list = mlr_var_list_int,
                                                   variable_list = variable_list_adj_log_iqr(),
                                                   conf_list = conf_list(),
                                                   sdq = "h_sdq_internal"),


  # We additionally adjusted our main model for breastfeeding and fish and seafood consumption
  # during pregnancy (since fish and seafood may accumulate persistent organic contaminants and
  # heavy metals).

  # For both analyses, the coefficients and CIs are similar as for the main ExWAS.

  # Breastfeeding adjustment
  # SDQ external
  bf_sens_ext = ExWAS::RunExwasImputed(input_data = imp_100_adj_log_iqr,
                                       variable_list = variable_list_ExWAS_ext,
                                       sdq = "h_sdq_external",
                                       conf_list = c(conf_list(), "breastf_cat"),
                                       save = TRUE,
                                       path = "results/supplementary_analyses",
                                       file_name = "sensitivity_breastf_ext"),

  # SDQ internal
  bf_sens_int = ExWAS::RunExwasImputed(input_data = imp_100_adj_log_iqr,
                                       variable_list = variable_list_ExWAS_int,
                                       sdq = "h_sdq_internal",
                                       conf_list = c(conf_list(), "breastf_cat"),
                                       save = TRUE,
                                       path = "results/supplementary_analyses",
                                       file_name = "sensitivity_breastf_int"),

  # Fish and seafood intake adjustment
  # SDQ external
  fish_sens_ext = ExWAS::RunExwasImputed(input_data = imp_100_adj_log_iqr,
                                         variable_list = variable_list_ExWAS_ext,
                                         sdq = "h_sdq_external",
                                         conf_list = c(conf_list(), "fish_intake"),
                                         save = TRUE,
                                         path = "results/supplementary_analyses",
                                         file_name = "sensitivity_fish_ext"),

  # SDQ internal
  fish_sens_int = ExWAS::RunExwasImputed(input_data = imp_100_adj_log_iqr,
                                         variable_list = variable_list_ExWAS_int,
                                         sdq = "h_sdq_internal",
                                         conf_list = c(conf_list(), "fish_intake"),
                                         save = TRUE,
                                         path = "results/supplementary_analyses",
                                         file_name = "sensitivity_fish_int"),

  # We explored sex-specific effects by adding an interaction term between each exposure and child sex.

  # Check if there is an interaction with sex with any of the exposures identified in the main ExWAS
  # (for one imputed dataset)

  # SDQ external
  sex_sens_ext = ExWAS::InteractionSex(input_data = imp_100_adj_log_iqr,
                                       sdq = "h_sdq_external",
                                       conf_list = conf_list(),
                                       strat_var_list = variable_list_ExWAS_ext$Variable_name) %>%

    # Display interactions with p value below 0.2
    dplyr::filter(p.value < 0.2), # none

  # SDQ internalizing
  sex_sens_int = ExWAS::InteractionSex(input_data = imp_100_adj_log_iqr,
                                       sdq = "h_sdq_internal",
                                       conf_list = conf_list(),
                                       strat_var_list = variable_list_ExWAS_ext$Variable_name) %>%

    # Display interactions with p value below 0.2
    dplyr::filter(p.value < 0.2), # none

  # We performed an ExWAS restricted to the participants with no missing biomarker concentrations.

  # IQR transform the non-imputed, standardized and log2 transformed dataset
  sdq_expo_adj_log_iqr = PrepareData::IqrTransform(dataset = sdq_expo_adj_log_cov),

  # SDQ external
  complete_case_ext = ExWAS::RunExwasCompleteCase(
    input_data = sdq_expo_adj_log_iqr,
    variable_list = variable_list_ExWAS_ext,
    sdq = "h_sdq_external",
    conf_list = conf_list()),

  # SDQ internal
  complete_case_int = ExWAS::RunExwasCompleteCase(
    input_data = sdq_expo_adj_log_iqr,
    variable_list = variable_list_ExWAS_int,
    sdq = "h_sdq_internal",
    conf_list = conf_list()),


  # For the exposures associated with the SDQ externalizing score we ran an ExWAS after exclusion
  # of the BiB cohort, as we had noted that children from this population had markedly lower
  # externalizing score (median = 0.5) compared to the other cohorts (medians ≥ 5).

  # Remove BiB cohort from the study
  BiB_dataset = dplyr::filter(stacked_dataset, cohort != "BIB") %>%
    droplevels(),

  # Calculate n for a new dataset withour BiB
  n_no_BiB = BiB_dataset %>%
    dplyr::filter(.imp == 0) %>%
    nrow(),

  # Collapse the new dataset into a mids object
  stacked_data_no_BiB = mice::as.mids(BiB_dataset),

  # SDQ external
  exwas_ext_no_BiB = ExWAS::RunExwasImputed(
    input_data = stacked_data_no_BiB,
    variable_list = variable_list_ExWAS_ext,
    sdq = "h_sdq_external",
    conf_list = conf_list()),

  # Because excessive maternal weight gain during pregnancy could lead to decreased blood
  #concentrations of lipophilic compounds due to their storage in the adipose tissue
  # [@kim_prenatal_2011; @lee_chlorinated_2014; @verner2013] and to behavioural problems
  # in the offspring [@pugh_gestational_2016], we ran an additional analysis stratified
  # on gestational weight gain for all the biomarkers from the OCs family.

  # Stratify ExWAS results on gestational weight gain in reference to pre-pregnancy BMI for compounds from the OC family
  OCs_stratified_weight_gain = ExWAS::StratifyGainWeight(data = complete_dataset,
                                                         variable_list = variable_list_adj_log_iqr(),
                                                         strat_var_list = variable_list_oc()),

  # Analyses included in the text only ----

  # Child birth dates

  # Display the birth dates of children
  child_age = table(sdq_expo_cov$year_birth, sdq_expo_cov$cohort) %>%
    print(),

  # Correlations between exposure concentrations

  # Note: For the study of the correlation between exposure concentrations, non-transformed, original
  # values of the concentrations will be used.

  # Prepare the list of exposures (remove Tl and DEDTP)
  # Select original untransformed values of exposures and change names for plotting
  sdq_expo_cov_corr = dplyr::select(sdq_expo_cov, ends_with("_m")) %>%
    data.table::setnames(old = colnames(.),
                         new = dplyr::filter(variable_list_m(), Exposure %nin% c("Tl", "DEDTP")) %>%
                           dplyr::pull(Exposure)),

  # Calculate correlations
  corr = DescriptiveStats::CalculateCorrelations(sdq_expo_cov_corr) %>%
    print(),

  # Geometrical mean for copper

  # Calculate geometric mean for copper
  geomean_Cu = DescriptiveStats::CalculateGeomean(sdq_expo_cov$hs_cu_m) %>%
    print(),


  # Paper Tables ----

  # Table 1: Population characteristics for the mother-child pairs included in the study: overall and by cohort
  res = read.csv(here::here("results/main_analyses", "overall_population_charactersitics.csv")) %>%
    as.data.frame(),

  res_coh = read.csv(here::here("results/main_analyses", "per_cohort_population_charactersitics.csv")) %>%
    as.data.frame() %>%
    dplyr::select(-X, ),

  Table_1 = Tab_1(res = res,
                  res_coh = res_coh,
                  path = "results/main_analyses",
                  file_name = "tables/Table_1"),

  # Table 2: Adjusted associations between the prenatal exposure to environmental contaminants
  # and SDQ externalising and internalising scores

  # Merge results from all analyses

  Table_2 = Tab_2(exwas_ext = exwas_ext,
                  exwas_int = exwas_int,
                  MLR_ext = MLR_ext,
                  MLR_int = MLR_int,
                  complete_case_ext = complete_case_ext,
                  complete_case_int = complete_case_int,
                  expo_sdq_ext = expo_sdq_ext,
                  expo_sdq_int = expo_sdq_int,
                  exwas_ext_no_BiB = exwas_ext_no_BiB,
                  n_no_BiB = n_no_BiB,
                  path = "results/main_analyses",
                  file_name = "tables/Table_2"),

  # Paper Figures ----

  # Figure 1 Cohort-stratified analysis of the associations between prenatal exposures
  # and SDQ external (A) or internal (B) scores run for each prenatal exposure variable
  # detected in the ExWAS (uncorrected p values <0.05), with an interaction term between
  # exposure and cohort origin added.

  # SDQ external
  Figure_1A = Heterogeneity::PlotHeterogeneity(data = stacked_dataset,
                                               outcome = "h_sdq_external",
                                               conf_list = conf_list(),
                                               variable_list = variable_list_ExWAS_ext,
                                               strat_variable = "cohort",
                                               width = 2500,
                                               height = 2000,
                                               res = 200,
                                               file_name = "figures/Figure_1A",
                                               path = "results/main_analyses"),

  # SDQ external
  Figure_1B = Heterogeneity::PlotHeterogeneity(data = stacked_dataset,
                                               outcome = "h_sdq_internal",
                                               conf_list = conf_list(),
                                               variable_list = variable_list_ExWAS_int,
                                               strat_variable = "cohort",
                                               width = 2500,
                                               height = 666,
                                               res = 135,
                                               file_name = "figures/Figure_1B1",
                                               path = "results/main_analyses"),


  # Figure 2: Sensitivity analysis for exposure-SDQ externalising score associations
  # stratified on gestational weight gain.

  # Plot the results of weight gain stratification for OCs family
  ExWAS::PlotStratifiedResults(data_for_plotting = OCs_stratified_weight_gain,
                               colors = c("black", "green", "blue"),
                               labels = c("Adequate", "Excessive", "Insufficient"),
                               ylim = c(0.5, 1.8),
                               legend_pos = c(0.5, 0.85),
                               path = "results/main_analyses",
                               file_name = "figures/Figure_2"),


  # Paper Appendix Tables ----

  # Appendix Table 5: Descriptive statistics of all prenatal exposures. Distribution indicators
  # are given for raw (i.e. non-imputed and non-transformed) exposure values.

  Appendix_Table_5 = App_Tab_5(exposure_char = exposure_char,
                               variable_list_m = variable_list_m(),
                               path = "results/supplementary_analyses",
                               file_name = "tables/Appendix_Table_5"),

  # Appendix Table 6: Adjusted associations between prenatal exposure to environmental contaminants
  # and SDQ externalising score (n = 708).

  Appendix_Table_6 = App_Tab_6_7(exwas = exwas_ext,
                                 heterog_estimates = heterog_estimates_ext,
                                 path = "results/supplementary_analyses",
                                 file_name = "tables/Appendix_Table_6"),

  # Appendix Table 7: Adjusted associations between prenatal exposure to environmental contaminants
  # and SDQ internalising score (n = 708).

  Appendix_Table_7 = App_Tab_6_7(exwas = exwas_int,
                                 heterog_estimates = heterog_estimates_int,
                                 path = "results/supplementary_analyses",
                                 file_name = "tables/Appendix_Table_7"),

  # Paper Appendix Figures ----


  # Appendix Figure 1: Study flowchart.

  Appendix_Figure_1 = App_Fig_1(path = "results/supplementary_analyses",
                                file_name = "figures/App_Figure_1"),

  # Appendix Figure 2: GAM models with restricted cubic splines function fitted on the
  # log2 and IQR transformed prenatal concentration of exposure selected in the LASSO
  # and ExWAS as associated with the SDQ score.

  Appendix_Figure_2A = App_Fig_2(GAM = GAM_ext,
                                 variable_list = variable_list_adj_log_iqr(),
                                 expo_sdq = expo_sdq_ext,
                                 path = "results/supplementary_analyses",
                                 file_name = "figures/App_Figure_2A"),

  Appendix_Figure_2B = App_Fig_2(GAM = GAM_int,
                                 variable_list = variable_list_adj_log_iqr(),
                                 expo_sdq = expo_sdq_int,
                                 path = "results/supplementary_analyses",
                                 file_name = "figures/App_Figure_2B"),

  # Answer to reviewer ----

  # ENET analysis

  # SDQ external
  enet_ext = LassoENET::ENET(input_data = stacked_dataset_imp_only,
                             sdq = "h_sdq_external",
                             variable_list = variable_list_adj_log_iqr(),
                             conf_list = conf_list(),
                             penalty_factor = penalty_factor(),
                             path = "results/revision",
                             file_name = "enet_ext"),

  # SDQ internal
  enet_int = LassoENET::ENET(input_data = stacked_dataset_imp_only,
                             sdq = "h_sdq_internal",
                             variable_list = variable_list_adj_log_iqr(),
                             conf_list = conf_list(),
                             penalty_factor = penalty_factor(),
                             path = "results/revision",
                             file_name = "enet_int"),

  # Display which variables exceeded the 50% selection threshold
  # SDQ external
  hits_enet_ext = LassoENET::HitsNumber(regr_result = enet_ext,
                                        variable_list = variable_list_adj_log_iqr()),

  # SDQ internal
  hits_enet_int = LassoENET::HitsNumber(regr_result = enet_int,
                                        variable_list = variable_list_adj_log_iqr()),


  # Revision Figure 1

  # SDQ external
  Revision_Figure_1A = LassoENET::PlotNoHits(no_hits = hits_enet_ext,
                                              title = "SDQ externalising",
                                              path = "results/revision",
                                              file_name = "enet_ext"),

  # SDQ internal
  Revision_Figure_1B = LassoENET::PlotNoHits(no_hits = hits_enet_int,
                                             title = "SDQ internalising",
                                             path = "results/revision",
                                             file_name = "enet_int"),

  # ExWAS on variables selected by ENET
  Revision_Table_1 = ExWAS::RunExwasSdqConfounders(input_data = imp_100_adj_log_iqr,
                                                    sdq = "h_sdq_external",
                                                    vars_to_regress = dplyr::pull(
                                                      dplyr::filter(
                                                        hits_enet_ext, number_of_hits >= 50),
                                                      Variable_name),
                                                    variable_names = variable_list_adj_log_iqr(),
                                                    conf_list = conf_list(),
                                                    path = "results/revision",
                                                    file_name = "Revision_Table_1"),

  # # * Main results ----
  # main_results = rmarkdown::render(
  #   drake::knitr_in("jedynak_prenatal_2020.Rmd"),
  #   output_file = drake::file_out("jedynak_prenatal_2020.html"),
  #   quiet = TRUE),

  # # * Answer to reviewer ----
  # reviewer_answer = rmarkdown::render(
  #   drake::knitr_in("jedynak_prenatal_2020_answer_to_reviewer.Rmd"),
  #   output_file = drake::file_out("jedynak_prenatal_2020_answer_to_reviewer.html"),
  #   quiet = TRUE),

  # * Export session info ----
  session_info = writeLines(utils::capture.output(utils::sessionInfo()), "sessionInfo.txt")
)

