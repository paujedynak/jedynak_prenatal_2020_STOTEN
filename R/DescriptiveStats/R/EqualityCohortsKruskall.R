#' Calculates equality of exposure concentrations between groups (cohorts in this case) using Kruskall-Walis test
#'
#' @param x A vector of exposure concentrations for each patient
#' @param coh A vector of cohort assigned to each patient
#'
#' @importFrom stats kruskal.test
#'
#' @return A scalar with a p value for the Kruskall-Walis test

.EqualityCohortsKruskall <- function(x, coh) {

  res <- stats::kruskal.test(x, coh)
  p_value <- res$p.value

  if (p_value < 0.0001) {
    p_value <- "< 0.0001"

  } else {
    p_value <- format(round(p_value, 3), nsmall = 3)
  }

  return(p_value)
}
