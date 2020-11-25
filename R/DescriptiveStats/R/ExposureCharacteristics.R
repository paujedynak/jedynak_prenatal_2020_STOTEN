# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("Variable_name",
                         "status",
                         "quant_LOD",
                         "below_LOD",
                         "above_LOD",
                         "perc_above_LOD",
                         "miss_nb",
                         "miss_perc",
                         "Concentration",
                         "Variable_short",
                         "cohort",
                         "median",
                         "iqr",
                         "cohort",
                         "Overall",
                         "V1"))


#' Calculates exposure characteristics
#'
#' @param input_data_descr Data frame with variables containing exposure characteristics (descriptive about LOD)
#' @param input_data_conc Data frame with variables containing exposure concentrations
#'
#' @return A data frame with exposure descriptive statistics
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom stats IQR
#' @importFrom stats median
#' @importFrom stringr str_replace
#' @importFrom tibble rownames_to_column
#'
ExposureCharacteristics <- function(input_data_descr, input_data_conc) {

# _Columns 3-4: Missing: Nb. (%) and  >LOD: Nb. (%)_

n <- nrow(input_data_descr)

miss_LOD <- input_data_descr %>%
  tidyr::gather(key = Variable_name, value = status) %>%
  dplyr::group_by(Variable_name) %>%

  # Count "Quantifiable value" and "Values <LOD" (all that are different than "Sample not analysed")
  # Count "Values <LOD"
  # Count "Values >LOD"
  # Count "perc_above_LOD"
  # Calculate percentage of "Values > LOD"
  dplyr::summarise(miss_nb = length(Variable_name[status == "Sample not analysed"]),
                   quant_LOD = length(Variable_name[status != "Sample not analysed"]),
                   below_LOD = length(Variable_name[status == "Values <LOD"]),
                   above_LOD = quant_LOD - below_LOD,
                   perc_above_LOD = above_LOD * 100 / quant_LOD,
                   '>LOD: Nb (%)' = paste0(above_LOD, " (", format(round(perc_above_LOD, 1), nsmall = 1), ")")) %>%
  dplyr::mutate(miss_perc = format(round(miss_nb * 100 / n, 1), nsmall = 1),
                "Missing: Nb. (%)" = paste0(miss_nb, " (", miss_perc, ")"),
                Variable_short = stringr::str_replace(Variable_name, "_mdesc", "")) %>%
  dplyr::select(Variable_short, 'Missing: Nb. (%)', '>LOD: Nb (%)')

# _Column 5: Median (IQR) (Overall)_

median_IQR <- input_data_conc %>%
  tidyr::gather(key = Variable_name, value = Concentration, -cohort) %>%
  dplyr::group_by(Variable_name) %>%
  dplyr::summarise(median = round(stats::median(Concentration, na.rm = TRUE), 1),
                   iqr = round(stats::IQR(Concentration, na.rm = TRUE), 1),
                   Overall = paste0(median, " (", iqr, ")")) %>%
  dplyr::mutate(Variable_short = stringr::str_replace(Variable_name, "_m$", "")) %>%
  dplyr::select(Variable_short, Overall)


# _Columns 6-10: Median (IQR) (per cohort)_

median_IQR_per_cohort <- input_data_conc %>%
  tidyr::gather(key = Variable_name, value = Concentration, -cohort) %>%
  dplyr::group_by(cohort, Variable_name) %>%
  dplyr::summarise(median = format(round(stats::median(Concentration, na.rm = TRUE), 1), nsmall = 1),
                   iqr = format(round(IQR(Concentration, na.rm = TRUE), 1), nsmall = 1),
                   'Median (IQR) cohort' = paste0(median, " (", iqr, ")")) %>%
  dplyr::mutate(Variable_short = stringr::str_replace(Variable_name, "_m$", "")) %>%
  dplyr::select(cohort, Variable_short, 'Median (IQR) cohort') %>%
  tidyr::spread(key = cohort, value = 'Median (IQR) cohort') %>%
  dplyr::mutate_all(~ stringr::str_replace(., stringr::fixed("NA (NA)"), "NA"))

# _Column 11: p value of equality between cohorts_

coh <- input_data_conc$cohort

data_equality <- input_data_conc %>%
  dplyr::summarise(dplyr::across(.cols = -cohort, .fns = ~ .EqualityCohortsKruskall(.x, coh = coh))) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Variable_name") %>%
  dplyr::rename("p value of equality between cohorts" = V1) %>%
  dplyr::mutate(Variable_short = stringr::str_replace(Variable_name, "_m$", "")) %>%
  dplyr::select(Variable_short, "p value of equality between cohorts")

# Create an object containing all the results
exposure_char <- list("miss_LOD" = miss_LOD,
                      "median_IQR_overall" = median_IQR,
                      "median_IQR_per_cohort" = median_IQR_per_cohort,
                      "p_value_Kruskall-Wallis" = data_equality)

return(exposure_char)
}
