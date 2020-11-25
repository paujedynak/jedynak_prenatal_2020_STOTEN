# Figures for Jedynak et al., 2020 paper
# Paulina Jedynak
# 23/09/20

# Appendix Figure 1

App_Fig_1 <- function(path, file_name) {
  main_x <- 0.31
  main_x_out <- 0.75

  gp <- grid::gpar(fill = "white")

  lwd <- 1

  lty_gp <- grid::gpar(lwd = lwd)

  lev_1_text <- "Entire cohorts, n = 31,472\n (BiB: 10,849, EDEN: 1,900, INMA: 2,063, KANC: 4,107, MoBa: 11,095, RHEA: 1,458)"
  lev_2_text <- "Women enrolled in the HELIX study, n = 1,301\n (BiB: 205, EDEN: 198, INMA: 223, KANC: 204, MoBa: 272, RHEA: 199)"
  lev_3_text <- "SDQ scores available, n = 721\n (BiB: 51, EDEN: 196, INMA: 221, KANC: 83, MoBa: 0, RHEA: 170)"
  lev_4_text <- "Study population, n = 708\n (BiB: 46, EDEN: 193, INMA: 218, KANC: 83, RHEA: 168)"

  out_3_text <- "Externalising and/or internalising SDQ scores not available, n = 580\n (BiB: 154, EDEN: 2, INMA: 2, KANC: 121, MoBa: 272, RHEA: 29)"
  out_4_text <- "Child sex or age at the SDQ assesment not available, n = 13\n (BiB: 5, EDEN: 3, INMA: 3, KANC: 0, RHEA: 2)"

  # TO VISUALIZE RUN IN THE CONSOLE, RUN IN THE CHUNK IT DOES NOT WORK!!

  grid::grid.newpage()

  grDevices::jpeg(here::here(path, paste0(file_name, ".jpg")), width = 3300, height = 1800, res = 300)

  # create boxes
  lev_1 <- Gmisc::boxGrob(lev_1_text, x = main_x, y = 0.9, box_gp = gp)
  lev_2 <- Gmisc::boxGrob(lev_2_text, x = main_x, y = 0.62, box_gp = gp)
  lev_3 <- Gmisc::boxGrob(lev_3_text, x = main_x, y = 0.34, box_gp = gp)
  lev_4 <- Gmisc::boxGrob(lev_4_text, x = main_x, y = 0.06, box_gp = gp)

  out_3 <- Gmisc::boxGrob(out_3_text, x = main_x_out, y = 0.48, box_gp = gp)
  out_4 <- Gmisc::boxGrob(out_4_text, x = 0.73, y = 0.2, box_gp = gp)

  Gmisc::connectGrob(lev_1, lev_2, "vertical", lty_gp = lty_gp)
  Gmisc::connectGrob(lev_2, lev_3, "vertical", lty_gp = lty_gp)
  Gmisc::connectGrob(lev_3, lev_4, "vertical", lty_gp = lty_gp)

  Gmisc::connectGrob(lev_3, out_3, "-", lty_gp = lty_gp)
  Gmisc::connectGrob(lev_4, out_4, "-", lty_gp = lty_gp)

  lev_1; lev_2; lev_3; lev_4; out_3; out_4

  grDevices::dev.off()

}

# Appendix Figure 2

App_Fig_2 <- function(GAM,
                      variable_list,
                      expo_sdq,
                      path,
                      file_name) {

  # Plot result for SDQ external
  grDevices::jpeg(here::here(path, paste0(file_name, ".jpg")), width = 2500, height = 2000, res = 300)

  plot_size <- length(GAM)
  graphics::par(mfrow = c(ceiling(plot_size / 2), 2))

  GAMrcs::PlotGAM(data = GAM,
                  variable_list = dplyr::filter(variable_list,
                                                Variable_name %in% expo_sdq))

  grDevices::dev.off()
}

