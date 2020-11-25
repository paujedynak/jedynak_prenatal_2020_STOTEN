# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("Exposure", "exp_lev"))


#' Calculate how many times the lasso selected each variable
#' @param regr_result A dataframe containing coefficients for each exposure for each run
#' @param variable_list A data.frame containing a list of exposures for which regressions will be run, first column must contain exposure names
#'
#' @return A dataframe with number of hits for each exposure variable
#'
#' @import dplyr
#' @importFrom magrittr set_colnames
#' @importFrom tibble rownames_to_column
#'
#' @export


HitsNumber <- function(regr_result, variable_list) {

  no_hits <- colSums(regr_result != 0) %>%
    data.frame() %>%
    magrittr::set_colnames("number_of_hits") %>%
    tibble::rownames_to_column("Variable_name") %>%
    merge(variable_list, ., by = "Variable_name") %>%
    dplyr::mutate(exp_lev = factor(Exposure, levels = Exposure),
                  Exposure = factor(exp_lev, levels = rev(levels(exp_lev))))

  no_hits_exc_50 <- dplyr::filter(no_hits, number_of_hits >= 50)

  if (nrow(no_hits_exc_50) > 0) {
    cat("Exposures that exceeded the 50% selection threshold:\n")

    for( i in seq_len(nrow(no_hits_exc_50))) {
      print(as.character(no_hits_exc_50[i, "Exposure"]))
    }
  } else {
    cat("Exposures that exceeded the 50% selection threshold: none")
  }

  return(no_hits)
}
