#' Calculates Spearman's correlation
#'
#' @param data A matrix of exposure concentrations
#'
#' @return A matrix of correlation coefficients
#' @export
#' @importFrom stats cor

CalculateCorrelations <- function(data) {

  corr <- data %>%
    stats::cor(use = "pairwise.complete.obs",
               method = "spearman") %>%
    round(2) %>%
    format(nsmall = 2)

  return(corr)
}
