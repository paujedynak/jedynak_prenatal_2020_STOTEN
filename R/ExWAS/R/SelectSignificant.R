# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables("Variable_name")

#' Selects variables that were significant on desired p value level
#'
#' @param data A dataframe which is a result of ExWAS analysis
#' @param pvalue A scalar defining p value threshold of significance
#'
#' @return A character vector containing names of significant exposures
#' @export

SelectSignificant <- function(data, pvalue) {
  sign <- data %>%
    dplyr::filter(p_value < pvalue) %>%
    dplyr::pull(Variable_name) %>%
    as.character()

  return(sign)
}
