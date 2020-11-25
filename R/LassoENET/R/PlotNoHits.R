# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("Exposure", "number_of_hits"))

#' Plot the number of the lasso hits per variable
#'
#' @param no_hits A dataframe with number of hits for each exposure variable
#' @param title A string defining the title of the plot
#' @param path A string defining the path for saving the file
#' @param file_name A string defining file name (without extension)
#'
#' @return A barplot of number of hits per variable
#'
#' @import ggplot2
#' @export

PlotNoHits <- function(no_hits, title, path, file_name) {

  plot_enet <- ggplot2::ggplot(no_hits, ggplot2::aes(Exposure, number_of_hits)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_minimal() +
    ggplot2::geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
    ggplot2::scale_y_continuous(name = "Number of hits (%)", limits = c(0, 100), breaks = seq(0, 100, 10), expand = c(0, 0)) +
    ggplot2::scale_x_discrete(name = "", expand = c(0, 0)) +
    ggplot2::coord_flip() +
    ggplot2::ggtitle(title)

  # save the plot to a file
  ggplot2::ggsave(here::here(path, paste0(file_name, ".jpg")),
                  plot = ggplot2::last_plot(),
                  device = "jpeg",
                  scale = 1,
                  width = 15,
                  height = 20,
                  units = "cm",
                  dpi = 150,
                  limitsize = TRUE)

  return(plot_enet)
}
