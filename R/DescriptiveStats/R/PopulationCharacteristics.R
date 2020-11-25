#' Calculates population characteristics
#'
#' @param input_data Data frame with variables containing population characteristics
#' @param variable_list A character vector with names of population characteristics variables
#' @param path A string defining the path to save the file
#'
#' @return A dataframe with population characteristics
#' @import compareGroups
#' @import dplyr
#' @importFrom here here
#' @importFrom tidyselect all_of
#' @export

PopulationCharacteristics <- function(input_data, variable_list, path) {

  # Statistics irrespectable of the cohort
  res <- compareGroups::compareGroups(data = dplyr::select(input_data,
                                                           tidyselect::all_of(variable_list), -cohort),
                                      include.miss = TRUE,
                                      method = c(h_sdq_external = 2,
                                                 h_sdq_internal = 2,
                                                 mother_age = 2,
                                                 child_age = 2)) %>%
    compareGroups::createTable(digits = (child_age = 1),
                               show.n = FALSE)

  # Statistics by cohort
  res_coh <- compareGroups::compareGroups(formula = cohort ~ .,
                                          data = dplyr::select(input_data,
                                                               tidyselect::all_of(variable_list)),
                                          include.miss = TRUE,
                                          method = c(mother_age = 2,
                                                     child_age = 2,
                                                     h_sdq_external = 2,
                                                     h_sdq_internal = 2)) %>%
    compareGroups::createTable(digits = (child_age = 1))

  # Save to file
  compareGroups::export2csv(res,
                            file = here::here(path, "overall_population_charactersitics.csv"),
                            header.labels = c(p.overall = "p value"))

  compareGroups::export2csv(res_coh,
                            file = here::here(path, "per_cohort_population_charactersitics.csv"),
                            header.labels = c(p.overall = "p value"))

  result <- list(Overall = res,
                 Per_cohort = res_coh)

  return(result)
}
