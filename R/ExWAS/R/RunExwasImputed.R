#' Run an ExWAS on an inputed `mids` object
#'
#' @param input_data "mids" object (imputed dataset)
#' @param variable_list A data.frame containing a list of exposures for which regressions will be run, first column must contain exposure names
#' @param sdq A character vector indicating in the outcome SDQ score is external or internal
#' @param subcohort A character vector specifying stratification to be made (default = NULL)
#' @param conf_list A character vector specifying confounders to be added in the regression model
#' @param save A character vector specifying if the files are going to be saved (default = FALSE, do not save files)
#' @param path Path for saving (default = NULL)
#' @param file_name File name (to which .csv or _plot.csv will be added, default = NULL)
#'
#'
#' @return Regression coefficients (IRR) with 95% confidence intervals for imputed dataset
#' @export
#'
RunExwasImputed <- function(input_data,
                            variable_list,
                            sdq,
                            conf_list,
                            subcohort = NULL,
                            save = FALSE,
                            path = NULL,
                            file_name = NULL) {

  # Run regressions
  exwas_result <- .RunRegression(input_data,
                                 variable_list,
                                 sdq,
                                 subcohort,
                                 conf_list)

  # Tidy results for all exposures
  exwas_tidy_result <- .TidyExwasAllOutcomes(input_data,
                                             variable_list,
                                             exwas_result,
                                             subcohort,
                                             save,
                                             path,
                                             file_name)
  return(exwas_tidy_result)
}
