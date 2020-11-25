# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables("Description")

#' Plot GAMs for a list of exposure
#'
#' @param data A list of GAM model(s)
#' @param variable_list A data.frame containing a list of exposures for wich regressions will be run, first column must contain exposure names
#' @param residuals A logical indicating if to display residuals or not (default = FALSE)
#'
#' @return A panel of plots
#' @export
#' @import dplyr
#'
PlotGAM <- function(data, variable_list, residuals = FALSE) {

  for (i in seq_along(data)) {
    plot_title <- variable_list$Description[i]

    plot <- plot(data[[i]],
                 se = TRUE,
                 terms = paste0("rms::rcs(", names(data)[i], ", knots = 3)"),
                 pch = 1,
                 cex = 0.75,
                 xlab = "Exposure biomarker concentration",
                 ylab = "Mean change in SDQ score",
                 main = plot_title,
                 ylim = c(-2, 2),
                 residuals = residuals)
  }

  return(plot)
}
