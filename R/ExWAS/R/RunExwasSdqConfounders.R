# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("term", "Variable_name", "Description"))

#' Runs ExWAS for one SDQ score at a time and a confounder
#'
#' @param input_data "mids" object (imputed dataset)
#' @param sdq A character vector indicating in the outcome SDQ score is external or internal
#' @param conf_list A character vector containing names of adjusting factors
#' @param vars_to_regress character vector containing names of variables that will be included in the multiple regression
#' @param variable_names A data.frame containing a list of exposures for which regressions will be run, first column must contain exposure names
#' @param file_name Name under which the result will be saved (*.csv)
#' @param path Path for saving
#'
#' @return Regression coefficients (IRR) with 95% confidence intervals for one confounder vs. SDQ score at a time
#' @export
#' @import mice
#' @import dplyr
#' @import utils
#' @importFrom MASS glm.nb
#' @importFrom stringr str_c
#' @importFrom tibble rownames_to_column
#' @importFrom here here
#' @importFrom tidyselect all_of
#'

RunExwasSdqConfounders <- function(input_data,
                                   sdq,
                                   vars_to_regress = NULL,
                                   conf_list,
                                   variable_names,
                                   path,
                                   file_name) {

  tidy_result <- list()

  # Change the 'mids' object imputed dataset to complete, stacked dataset
  all_imps <- mice::complete(input_data, "long") %>%

    # Select variables of interest plus covariates, SDQ cores and imp number
    dplyr::select(.imp, tidyselect::all_of(sdq), tidyselect::all_of(vars_to_regress), tidyselect::all_of(conf_list))

  if (!is.null(vars_to_regress)) {

    formula <- stats::as.formula(paste(sdq, "~",
                                       paste(vars_to_regress, collapse = " + "),
                                       " + ",
                                       paste(conf_list, collapse = " + ")))

    # perform regression for each imputed dataset (group by
    # imputation number)
    pool_fit <- all_imps %>%

      dplyr::group_by(.imp) %>%

      # run simple regression for sdq ~ exposure + covariates
      dplyr::do(model = MASS::glm.nb(formula = formula, data = .))

  } else {

    for (conf in conf_list) {

      #If confounder is a cohort, fit the model with SDQ and cohort only
      if (conf == "cohort") {

        # change input dataset to dataframe
        pool_fit <- all_imps %>%

          # perform regression for each imputed dataset (group by
          # imputation number)
          dplyr::group_by(.imp) %>%

          # run simple regression for sdq ~ exposure + covariates
          dplyr::do(model = MASS::glm.nb(formula = get(sdq) ~ get(conf), data = .))

      } else {

        # If confounder is not a cohort, calculate how many levels the confounder has
        nlev <- dplyr::select(all_imps, conf) %>%
          .[, 1] %>%
          levels() %>%
          length()


        # fit the model with SDQ, cohort and confounder
        pool_fit <- all_imps %>%
          dplyr::group_by(.imp) %>%
          dplyr::do(model = MASS::glm.nb(formula = get(sdq) ~ get(conf) + cohort, data = .))
      }
    }
  }

  # transform the regression output to a list
  tidy_fit <- pool_fit %>%
    as.list() %>%

    # pick the list item containing numerical regression output
    # (coefficient and other statistics for the intercept,
    # confounders and exposure), ignore other results provided by
    # the glm.nb function
    .[[-1]] %>%

    # pool the regression outcomes for each imputed dataset, create new columns containing estimates
    #and CIs, assign name to the row with the confounder
    mice::pool() %>%
    summary(conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::rename(Variable_name = term,
                  conf_low = '2.5 %',
                  conf_high = '97.5 %',
                  p_value = p.value)

  if (!is.null(vars_to_regress)) {

    tidy_fit <- dplyr::filter(tidy_fit, Variable_name %in% vars_to_regress)

  } else {

    if (conf == "cohort") {
      tidy_fit <- tidy_fit[-1, ]

    } else if (class(dplyr::select(all_imps, conf)[, 1]) == "factor") {

      tidy_fit <- tidy_fit[2:nlev, ]

    } else {

      tidy_fit <- tidy_fit[2, ]
    }
  }

  tidy_fit <- tidy_fit %>%
    dplyr::mutate(CI = stringr::str_c(round(conf_low, 2), round(conf_high, 2), sep = "; "),
                  estCI = stringr::str_c(round(estimate, 2), " (", CI, ")", sep = ""),
                  p_value = round(p_value, 3)) %>%
    dplyr::select(-CI)

  if (is.null(vars_to_regress)) {

    tidy_fit <- tidy_fit %>%
      dplyr::mutate(coadj = conf) %>%
      dplyr::mutate(confounder = gsub("get[(]conf[)]", paste0(conf, "_"), confounder)) %>%
      dplyr::select(confounder, coadj, estCI, p_value)

    tidy_result[[length(tidy_result) + 1]] <- tidy_fit

  } else {

    tidy_result[[length(tidy_result) + 1]] <- tidy_fit
  }

  # Combine results into one dataset
  res_all <- as.data.frame(do.call(rbind, tidy_result)) %>%
    merge(variable_names, . , by = "Variable_name") %>%
    dplyr::select(Description, Family, estCI, p_value)

  # Save result to a file
  utils::write.csv(res_all, here::here(path, paste0(file_name, ".csv")))

  return(res_all)
}
