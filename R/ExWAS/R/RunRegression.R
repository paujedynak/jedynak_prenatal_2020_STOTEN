# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("sex", "cohort"))

#' Run adjusted regression
#'
#' @param input_data "mids" object (imputed dataset)
#' @param variable_list A data.frame containing a list of exposures for which regressions will be run, first column must contain exposure names
#' @param sdq A character vector indicating in the outcome SDQ score is external or internal
#' @param subcohort A character vector specifying stratification to be made (default = NULL)
#' @param conf_list A character vector specifying confounders to be added in the regression model
#' @return A data.frame with regression results for all exposures in a non-tidy format
#'
#' @import dplyr
#' @import utils
#' @importFrom Helpers Mids2Stacked
#'


.RunRegression <- function(input_data,
                           variable_list,
                           sdq,
                           subcohort = NULL,
                           conf_list) {

  exwas_result <- data.frame()

  # Obtain stacked complete imputed dataset
  all_imps <- Helpers::Mids2Stacked(input_data, variable_list)

  if (is.null(subcohort) == FALSE) {

    # Extract subcohorts names
    subcohort_level <- levels(dplyr::pull(all_imps, subcohort))

    for (element in subcohort_level) {

      if (subcohort == "sex") {

        # Filter datasets leaving one sex at a time
        dataset <- dplyr::filter(all_imps, sex == element)

      } else if (subcohort == "cohort") {

        # Filter datasets removing one cohort at a time
        dataset <- dplyr::filter(all_imps, cohort != element)
      }

      # Remove unnecessary levels
      dataset <- droplevels(dataset)

      # Obtain pooled fits for regression fits all exposures
      tidy_res <- .ObtainPooledFits(dataset = dataset, variable_list = variable_list, sdq = sdq, conf_list = conf_list)

      if (subcohort == "sex") {

        tidy_res <- tidy_res %>%
          dplyr::mutate(sub_cohort = element)

      } else if (subcohort == "cohort") {

        tidy_res <- tidy_res %>%
          dplyr::mutate(sub_cohort = paste0("No_", element))
      }

      # Save results for each exposure
      exwas_result <- rbind(exwas_result, tidy_res)
    }

  } else {

    # Obtain pooled fits for regression fits all exposures
    tidy_res <- .ObtainPooledFits(dataset = all_imps, variable_list = variable_list, sdq = sdq, conf_list = conf_list)

    # Save results for each exposure
    exwas_result <- rbind(exwas_result, tidy_res)
  }

  return(exwas_result)
}

