# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

#' Impute missing values using mice package
#'
#' @param dataset A dataframe containing missing values
#' @param no_pred_list A character vector containing names of the variables to be skipped from imputation
#' @param high_corr_list A character vector containing names of the highly correlated variables to be treated specially
#' @param low_pred_list A character vector containing names of the variables with low prediction potential, to be treated specially
#' @param skip_list A character vector containing names of the variables to be skipped from imputation but not from prediction
#' @param creat_list A character vector containing names of the variables that were creatinine adjusted (passive imputation)
#' @param fat_list A character vector containing names of the variables that were total fat percentage adjusted (passive imputation)
#' @param log2_list List variables not adjusted for creatinine or fat but log2 transformed
#' @param m Number of imputations
#' @param seed Seed value
#' @param path A string with a path to save the file
#' @param file_name A string with a file name
#'
#' @return A mids object
#'
#' @export
#'
#' @import mice
#' @import dplyr
#' @importFrom here here

MultipleImputation <- function(dataset,
                               no_pred_list,
                               high_corr_list,
                               low_pred_list,
                               skip_list,
                               creat_list,
                               fat_list,
                               log2_list,
                               m,
                               seed,
                               path,
                               file_name) {


  # A few operations must be performed on the dataset in order to recreate the dataset
  # imputed 16.09.2019 that was used in the paper. Imputation is very sensitive for even
  # coding of variables, so to replicate the imputed dataset, the input datasets must be identical.
  # These operations will be reversed after the imputation.

  dataset_edit <- dataset %>%
    dplyr::mutate(cohort = factor(cohort,
                                  levels = c("EDEN", "BiB", "KANC", "RHEA", "INMA"),
                                  labels = c("EDEN", "BIB", "KANC", "RHEA", "INMA")),
                  mother_edu = factor(mother_edu,
                                      levels = c("Primary", "Secondary", ">=University"),
                                      labels = c("Low", "Middle", "High")),
                  mother_work = factor(mother_work, levels = c("Unemployed", "Employed"),
                                       labels = c("No", "Yes")))

  # 1. print initial predictor matrix

  init <- mice::mice(dataset_edit, maxit = 0)
  meth <- init$method
  predM <- init$predictorMatrix


  # 2. Change the predictor matrix

  # Remove HelixID and birth year from being a predictor
  predM[, no_pred_list] <- 0

  # Remove highly correlated variables as predictors for all the variables
  predM[, high_corr_list] <- 0

  # Remove variables with low predictive ability (basing on outflux/influx plots)
  predM[, low_pred_list] <- 0




  # 3. Change the method matrix

  # skip variables from imputation but not from prediction
  meth[skip_list] <- ""

  ## ------------------------------------------------------------------------
  # Create dehp sum variables
  meth["dehp_sum_m"] <- "~ I((hs_mehp_m / 278) + (hs_mehhp_m / 294) + (hs_meohp_m / 292) + (hs_mecpp_m / 308))"

  # Remove sum variable from prediction of the summed variables
  predM[c("hs_mehp_m", "hs_mehhp_m", "hs_meohp_m", "hs_mecpp_m"), "dehp_sum_m"] <- 0

  meth["weight_gain_cat"] <- "~ I(dplyr::case_when(mother_bmi == 'Underweight' & weight_gain < 12.5 | mother_bmi == 'Normal_weight' & weight_gain < 11.5 | mother_bmi == 'Overweight' & weight_gain < 7 | mother_bmi == 'Obesity' & weight_gain < 5 ~ 'Insufficient', mother_bmi == 'Underweight' & weight_gain > 18 | mother_bmi == 'Normal_weight' & weight_gain > 16 | mother_bmi == 'Overweight' & weight_gain > 11.5 | mother_bmi == 'Obesity' & weight_gain > 9 ~ 'Excessive', TRUE & !is.na(weight_gain) & !is.na(mother_bmi) ~ 'Adequate'))"

  predM[c("weight_gain", "mother_bmi"), "weight_gain_cat"] <- 0


  # Create variables adjusted for creatinine (passive imputation). Follow formula for Helix biomarkers protocol p. 15

  for (item in creat_list) {
    meth[item[1]] <- item[2]
  }

  # Remove adjusted variables from prediction of creatinine and the unadjusted variables
  for (item in creat_list) {
    predM[c(item[3], "creatinine", item[4]), item[1]] <- 0
  }

  ####################


  # Create variables adjusted for fat (passive imputation). Follow formula for Helix biomarkers protocol p. 14

  for (item in fat_list) {
    meth[item[1]] <- item[2]
  }


  # Remove adjusted variables from prediction of fat and the unadjusted variables
  for (item in fat_list) {
    predM[c(item[3], "fat", item[4]), item[1]] <- 0
  }

  ####################


  # Create not adjusted variables log2 transformed

  for (item in log2_list) {
    meth[item[1]] <- item[2]
  }

  # Remove log2 variables from prediction of the untransformed variables
  for (item in log2_list) {
    predM[item[3], item[1]] <- 0
  }

  # 4. Run the imputation

  # 100 imputations with a seed set
  imp_100_adj_log_temp <- mice::mice(data = dataset_edit,
                                     method = meth,
                                     predictorMatrix = predM,
                                     m = m,
                                     seed = seed)

  # Reverse the coding that was changed at the beginning of this script
  imp_100_adj_log <- mice::complete(data = imp_100_adj_log_temp,
                                    action = "long",
                                    include = TRUE) %>%
    dplyr::mutate(cohort = factor(cohort,
                                  levels = c("EDEN", "BIB", "KANC", "RHEA", "INMA"),
                                  labels = c("EDEN", "BiB", "INMA", "KANC", "RHEA")),
                  mother_edu = factor(mother_edu,
                                      levels = c("Low", "Middle", "High"),
                                      labels = c("Primary", "Secondary", ">=University")),
                  mother_work = factor(mother_work, levels = c("No", "Yes"),
                                       labels = c("Unemployed", "Employed"))) %>%
    mice::as.mids()

  # For easier comparison with other studies, after imputation the exposure concentrations will be IQR transformed.
  # Apply IQR transformation
  imp_100_adj_log_iqr <- PrepareData::IqrTransform(dataset = imp_100_adj_log)

  # Create an object containing both IQR transformed and untransformed imputed datasets

  imputed_datasets <- list("imputed_data" = imp_100_adj_log,
                           "imputed_IQR_data" = imp_100_adj_log_iqr)

  # Save object to the file
  saveRDS(imputed_datasets, here::here(path, file_name))

  # NAs are introduced when the stacked dataset is coerced to mids using as.mids (and this is ok)
  return(imputed_datasets)
}

