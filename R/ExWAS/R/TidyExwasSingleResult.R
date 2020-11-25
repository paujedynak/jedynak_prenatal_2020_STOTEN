# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

#' Generate tidy ExWAS result from a regression on a single exposure
#'
#' @param x Pooled result of a regression run on a mids object
#'
#' @return Statistics (regression coefficient, 95% CI, p-value...) for the modeled variable of interest
#' @importFrom tibble rownames_to_column
#'


.TidyExwasSingleResult <- function(x) {

  # transform rownames into 1st row
  tidy_fit <- tibble::rownames_to_column(x, var = "Exposure") %>%

    # pick the row for an exposure only, ignore coef. for
    # intercept and confounders
    .[2, ]

  return(tidy_fit)
}
