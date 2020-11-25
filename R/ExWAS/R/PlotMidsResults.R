# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables("exp_lev")

#' Plots regressions estimates obtained in ExWAS
#'
#' @param data_for_plotting Results of the ExWAS analysis containing: Exposure, Estimate, conf_low and conf_high
#' @param title Title of the plot
#' @param plot_name Name that the plot will be saved under (of the format: "plot.jpg")
#' @param path Path where to save the plot
#' @param ... other arguments passed to ggplot argument
#'
#' @return A plot of the regression coefficients IRR with 95% confidence intervals
#' @export
#' @import dplyr
#' @import utils
#' @import ggplot2
#' @importFrom here here
#'


PlotMidsResults <-
  function(data_for_plotting, title, plot_name, path, ...) {

    # Transform exposure into factor (to be able to order plot axis by exosure Family name)
    data <- data_for_plotting %>%
      dplyr::mutate(exp_lev = factor(Exposure, levels = Exposure),
             Exposure = factor(exp_lev, levels = rev(levels(exp_lev))))

    # Change order of the variables to add SumDEHP after its compounds
    #data_for_plotting <- data_for_plotting[c(1:47, 53, 48:52), ]

    # Plot data
    plot <-
      ggplot2::ggplot(data, aes(Estimate, Exposure, color = Family, fill = Family)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
      geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.1, size = 1.1) +
      geom_point(size = 3, shape = 21, color = "black") +
      xlab("IRR") +
      ylab("") +
      theme_bw() +
      ggtitle(title) +
      theme(text = element_text(size = 16), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = "bold"), legend.position = "none")

    # save the plot to a file
    ggplot2::ggsave(here::here(path, plot_name), plot = ggplot2::last_plot(), device = "jpeg", scale = 1, width = 15, height = 20, units = "cm", dpi = 150, limitsize = TRUE)

    return(plot)
  }
