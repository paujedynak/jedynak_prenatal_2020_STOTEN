# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("."))
}

#' log transforms exposure variables concentrations
#'
#' @param data A dataframe containing variables to be log transformed
#' @param list_log2 A character vector with names of variables to be log2 transformed
#' @param list_ln A character vector with names of variables to be ln transformed
#'
#' @return A dataframe containing original and transformed variables
#' @import dplyr

.LogTransform <- function(data, list_log2, list_ln) {

    if (class(data) != "data.frame") {
        stop("data must be a data.frame!")
    }

    if (class(list_log2) != "character" | class(list_ln) != "character") {
        stop("list_log2 and list_ln must be a character vector")
    }

    # Log2 variables that were previously creatinine/ fat adjusted, except cotinine
    log2_adj <- dplyr::select(data, ends_with("_adj"), -tidyselect::all_of(list_ln)) %>%
        sapply(log2) %>%
        as.data.frame()
    colnames(log2_adj) <- paste0(colnames(log2_adj), "_log2")

    # Log2 variables that were not creatinine/ fat adjusted previously
    log2_nonadj <- dplyr::select(data, tidyselect::all_of(list_log2)) %>%
        sapply(log2) %>%
        as.data.frame()
    colnames(log2_nonadj) <- paste0(colnames(log2_nonadj), "_log2")

    # ln variables (cotinine)
    ln_var <- dplyr::select(data, tidyselect::all_of(list_ln)) %>%
        sapply(log) %>%
        as.data.frame()
    colnames(ln_var) <- paste0(colnames(ln_var), "_ln")

    # Merge log2 and ln transformed variables with the rest of the dataset
    adj_sdq_expo_log <- cbind(data, log2_adj, log2_nonadj, ln_var)

    return(adj_sdq_expo_log)
}
