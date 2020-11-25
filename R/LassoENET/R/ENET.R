#' ENET with negative binomial outcome repeated for each imputed dataset
#'
#' @param input_data A dataframe containing a complete stacked dataset
#' @param sdq A string defining the output variable ("h_sdq_external" or "h_sdq_internal")
#' @param variable_list A character vector containing a list of exposures to be included in the lasso
#' @param penalty_factor An integer vector with 1 meaning variables to be shrunk and 0 variables that are forced not to be shrunk (the number is a sum of variables and contrasts for categorical variables)
#' @param path A string defining the path for saving the file
#' @param file_name A string defining file name (without extension)
#' @param conf_list A character vector listing confounders to be used in the model
#'
#' @return A dataframe containing number of ENET "hits" for each variable
#' @export
#' @import dplyr
#' @import mpath
#' @importFrom here here
#' @importFrom stats as.formula
#' @importFrom utils write.csv
#' @importFrom stats coef
#'
ENET <- function(input_data,
                 sdq,
                 variable_list,
                 conf_list,
                 penalty_factor,
                 path,
                 file_name) {

  # List exposure variable names
  exposures <- variable_list$Variable_name

  # Select variables with appropriate SDQ as the outcome variable and transform into formula
  ENET_formula <- stats::as.formula(paste(sdq, "~",
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

    # Provide a sequence of alphas for cross-validation
    a <- c(c(1 %o% 10^(-10:-1)), seq(from = 0.2, to = 1, by = 0.1))

    search_res <- data.frame()

    for (item in a) {

      cv <- mpath::cv.glmregNB(formula = ENET_formula,
                               data = imp_subset,
                               nfolds = 10,
                               alpha = item,
                               plot.it = FALSE,
                               penalty.factor = penalty_factor,
                               parallel = TRUE,
                               n.cores = parallel::detectCores())

      # Maximized log-likelihood for prediciton
      cv_optim <- max(cv$cv)

      # Calculate lambda.1se (basing on the glinternet.cv function)
      bestIndex1Std <- which(cv$cv >= cv_optim - cv$cv.error[which.max(cv$cv)])

      # I want to choose lambda higher than optimal one (increase of lambda = less variables selected) but not worsening the cv error more than for a value of 1 standard error. That is why I choose the first value of lambda for error + 1SE (as lambdas are tested from the highest to the smallest): bestIndex1Std[1].
      if (length(bestIndex1Std) == length(cv$lambda)) {
        lambdaHat1Std <- lambdaHat <- cv$lambda.optim

      } else {
        lambdaHat1Std <- cv$lambda[bestIndex1Std[1]]
        lambdaHat <- cv$lambda.optim
      }

      search <- data.frame(lambdaHat = lambdaHat,
                           cv_optim = cv_optim,
                           lambdaHat1Std = lambdaHat1Std,
                           alpha = item)

      search_res <- rbind(search_res, search)
    }

    alpha_optim <- search_res[search_res$cv_optim == max(search_res$cv_optim), "alpha"]

    lambda_1SE <- search_res[search_res$cv_optim == max(search_res$cv_optim), "lambdaHat1Std"]

    ENET_imp <- mpath::glmregNB(formula = ENET_formula,
                                data = imp_subset,
                                lambda = lambda_1SE,
                                alpha = alpha_optim,
                                penalty.factor = penalty_factor)

    # Remove intercept
    coef_enet_1SE <- stats::coef(ENET_imp)[-1]

    coefs_100[[i]] <- coef_enet_1SE
  }

  # Pool all coefficients from 100 models into a data frame
  pooled_coefs_100 <- as.data.frame(do.call(rbind, coefs_100))

  # Save th result to a file
  saveRDS(pooled_coefs_100, file = here::here(path, paste0(file_name, ".RDS")))

  return(pooled_coefs_100)
}
