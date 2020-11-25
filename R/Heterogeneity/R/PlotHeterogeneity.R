# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("est", "se"))

#' Draw forest heterogeneity plots - CODE and COMMENTS come from Xavier Basagaña
#'
#' @param data A complete stacked dataset (mice object), including .imp = 0
#' @param outcome A string defining an outcome ("h_sdq_external" or "h_sdq_internal")
#' @param conf_list A vector containing all confounder variable names
#' @param variable_list A data.frame containing a list of exposures for wich regressions will be run, first column must contain exposure names
#' @param strat_variable A string describing stratifying variable (normally "cohort")
#' @param width A scalar defining width of the plot
#' @param height A scalar defining height of the plot
#' @param res A scalar defining resolution of the plot
#' @param path A string defining the path for saving the plot
#' @param file_name A string defining the file name for saving the plot
#'
#' @return A forest plot
#' @export
#' @import metaplus
#' @import graphics
#' @import grDevices
#' @importFrom here here

PlotHeterogeneity <- function(data,
                              outcome,
                              conf_list,
                              variable_list,
                              strat_variable,
                              width,
                              height,
                              res,
                              path,
                              file_name) {

# List names of exposure variables
selected_exposures <- variable_list$Variable_name

# List descriptions of exposure variables
exposure_names <- variable_list$Description

# The following function assumes the following order of the cohorts: "BIB"  "EDEN" "INMA" "KANC" "RHEA"
estimates_plot <- .CalculateI(data,
                              outcome,
                              conf_list,
                              selected_exposures,
                              strat_variable)


grDevices::jpeg(here::here(path, paste0(file_name, ".jpg")),
                width = width,
                height = height,
                res = res)

graphics::par(mfrow = c(ceiling(nrow(variable_list) / 2), 2))

# Plot forest plot for all exposures
final_plot <- .FinalForestPlot(selected_exposures = selected_exposures,
                               exposure_names = exposure_names,
                               outcome = outcome,
                               RES = estimates_plot$RES,
                               RESv = estimates_plot$RESv)
grDevices::dev.off()

return(final_plot)
}
