#' Helper function to automatically change several variables' names
#'
#' @param x A name of the variable whose name is to be changed
#' @param hs_m A logical indicating which name transformation to perform
#' @param m_log2 A logical indicating which name transformation to perform
#' @param adj_log2 A logical indicating which name transformation to perform
#'
#' @export
#' @return A string with changed variable name

ChangeVarsNames <- function(x, hs_m = FALSE, m_log2 = FALSE, adj_log2 = FALSE) {

  # Make functions to change variable names

  if (isTRUE(hs_m)) {

    fun <- function(x) {
      x <- paste0("hs_", x, "_m")
      return(x)
    }

    res <- fun(x)
  }

  if (isTRUE(m_log2)) {

    fun <- function(x) {
      x <- paste0("hs_", x, "_m_log2")
      return(x)
    }

    res <- fun(x)
  }

  if (isTRUE(adj_log2)) {

    fun <- function(x) {
      x <- paste0("hs_", x, "_adj_log2")
      return(x)
    }

    res <- fun(x)
  }

  return(res)
}
