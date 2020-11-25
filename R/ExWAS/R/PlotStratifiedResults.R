# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables("group")

#' Plots regressions estimates obtained in stratified ExWAS
#'
#' @param data_for_plotting Results of the ExWAS analysis containing: Exposure, Estimate, conf_low and conf_high
#' @param colors Provide colors to be used for each stratum
#' @param ylim A numeric vector defining y axis limits (c(min, max))
#' @param legend_pos A numeric vector defining legend position (c(x, y))
#' @param labels Provide names to be used for each stratum
#' @param path A string defining the path for saving the plot
#' @param file_name A string defining the file name for saving the plot
#'
#' @return A plot of stratified regression
#' @export
#' @import ggplot2
#' @import utils
#'


PlotStratifiedResults <-
  function(data_for_plotting, colors, labels, ylim, legend_pos, path, file_name) {

    plot <-

      # Create a main plot
      ggplot2::ggplot(
        data_for_plotting,
        ggplot2::aes(Exposure, Estimate)) +

      # Add line marking IRR = 1
      ggplot2::geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "gray") +

      # Mark IRR estimates as points with 95% CI
      ggplot2::geom_pointrange(
        ggplot2::aes(ymin = conf_low,
                     ymax = conf_high,
                     color = group),
        position = ggplot2::position_dodge(0.5),
        size = 0.4) +

      # Add labels
      ggplot2::labs(
        title = "",
        x = "",
        y = "Incidence rate ratio (IRR)",
        color = "Gestational weight gain:") +

      # Add colors and labels to the legend
      ggplot2::scale_color_manual(
        values = colors,
        labels = labels) +

      # Change ylim
      ggplot2::ylim(ylim) +

      ggplot2::theme_classic() +

      # Add legend
      ggplot2::theme(
        #axis.text.x = element_blank(),

        legend.position = legend_pos,
        panel.background = ggplot2::element_rect(
          colour = "black",
          size = 1),

        # legend.box.background = element_rect(
        #   color = "black",
        #   fill = "grey90",
        #   size = 1,
        #   linetype = "solid"),

        legend.key.size = ggplot2::unit(0.25, "cm"),
        legend.key.width = ggplot2::unit(1, "cm")) +

      # Save the plot to a file
      ggplot2::ggsave(filename = here::here(path, paste0(file_name, ".jpg")),
                      plot = ggplot2::last_plot(),
                      device = "jpeg",
                      scale = 1,
                      width = 14,
                      height = 11,
                      units = "cm",
                      dpi = 150,
                      limitsize = TRUE)

    return(plot)
  }
