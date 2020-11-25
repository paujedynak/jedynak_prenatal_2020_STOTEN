#' Lasso with negative binomial outcome repeated for each imputed dataset
#'
#' @param input_data A dataframe containing a complete stacked dataset
#' @param sdq A string defining the output variable ("h_sdq_external" or "h_sdq_internal")
#' @param variable_list A character vector containing a list of exposures to be included in the lasso
#' @param penalty_factor An integer vector with 1 meaning variables to be shrunk and 0 variables that are forced not to be shrunk (the number is a sum of variables and contrasts for categorical variables)
#' @param path A path for saving the file
#' @param file_name File name (without extension)
#' @param conf_list A character vector listing confounders to be used in the model
#'
#' @return A dataframe containing coefficients for each variable for each run
#' @export
#' @import dplyr
#' @import mpath
#' @importFrom here here
#' @importFrom stats as.formula
#' @importFrom utils write.csv
#' @importFrom stats coef

Lasso <- function(input_data,
                  sdq,
                  variable_list,
                  conf_list,
                  penalty_factor,
                  path,
                  file_name) {

  # List exposure variable names
  exposures <- variable_list$Variable_name

  # Create regression formula
  lasso_formula <- stats::as.formula(paste(sdq, "~",
                                           paste(exposures, collapse = " + "),
                                           " + ",
                                           paste(conf_list, collapse = " + ")))

  # Create a list containing results of 100 regressions (on each imputed dataset)
  coefs_100 <- list()

  for (i in seq_len(max(input_data$.imp))) {

    # Counter of progress
    print(paste0("Imputed dataset no. ", i, "/", max(input_data$.imp)))

    # Subset each imputed dataset
    imp_subset <- dplyr::filter(input_data, input_data$.imp == i)

    # Run a cross-validation regression to define lambda.1se
    cv <- mpath::cv.glmregNB(formula = lasso_formula,
                             data = imp_subset,
                             lambda = NULL,
                             nFolds = 10,
                             alpha = 1,
                             plot.it = FALSE,
                             penalty.factor = penalty_factor,
                             parallel = TRUE,
                             n.cores = parallel::detectCores())

    # Calculate lambda.1se (basing on the glinternet.cv function)
    bestIndex1Std <- which(cv$cv >= max(cv$cv) - cv$cv.error[which.max(cv$cv)])

    if (length(bestIndex1Std) == length(cv$lambda)) {

      lambdaHat1Std <- lambdaHat <- cv$lambda.optim

    } else {

      lambdaHat1Std <- cv$lambda[bestIndex1Std[1]]
      lambdaHat <- cv$lambda.optim
    }

    # Run neg bin model with calculated lambda.1se
    lasso_imp  <- mpath::glmregNB(formula = lasso_formula,
                                  data = imp_subset,
                                  lambda = lambdaHat1Std,
                                  alpha = 1,
                                  plot.it = FALSE,
                                  penalty.factor = penalty_factor)

    # Save coefficients to a list
    coefs_lasso <- stats::coef(lasso_imp)[-1] # Remove intercept
    coefs_100[[i]] <- coefs_lasso
  }

  # Pool all coefficients from 100 models into a data frame
  pooled_coefs_100 <- as.data.frame(do.call(rbind, coefs_100))

  # Save th result to a file
  saveRDS(pooled_coefs_100, file = here::here(path, paste0(file_name, ".rds")))
  write.csv(pooled_coefs_100, file = here::here(path, paste0(file_name, ".csv")))

  return(pooled_coefs_100)
}
