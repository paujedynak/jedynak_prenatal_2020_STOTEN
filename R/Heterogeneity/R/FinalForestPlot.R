#' Plots forest plot for all exposures
#'
#' @param selected_exposures A character vector listing names of exposure variables
#' @param exposure_names A character vector listing names of exposure variables that will be displayed on the plots
#' @param outcome A string defining an outcome ("h_sdq_external" or "h_sdq_internal")
#' @param RES A list of exposure-strat_variable interaction
#' @param RESv A list
#'
#' @return Forest plot for each exposure
#' @export
#'

.FinalForestPlot <- function(selected_exposures,
                             exposure_names,
                             outcome,
                             RES,
                             RESv) {

  P <- ceiling(length(selected_exposures) / 6)

  j = 1

  for (p in seq_len(P)) {
    listVARt <- as.character(selected_exposures[(p * 6 - 5):min(p * 6, length(selected_exposures))])
    listVARnames <- exposure_names[(p * 6 - 5):min(p * 6, length(exposure_names))]

    for (k in seq_along(listVARt)) {
      print(listVARt[k])
      temp_plot <- RESv[[listVARt[k]]]

      # you can modify rownames(temp_plot) if you want (I have added symbols relatively to the proportion of missing values, i.e. if prop is the list of proportion of missing values for selected_exposures

      # wh <- which(prop<0.8 & prop>=0.30)
      # if(length(wh)>0 )rownames(temp)[wh] <- paste(rownames(temp)[wh],"*")
      # wh <- which(prop<0.3 & prop>=0.10)
      # if(length(wh)>0 )rownames(temp)[wh] <- paste(rownames(temp)[wh],"**")
      # wh <- which(prop<0.10 )
      # if(length(wh)>0 )rownames(temp)[wh] <- paste(rownames(temp)[wh],"***")

      mag_meta_plot <- metaplus::metaplus(yi = exp(est), sei = se, slab = rownames(temp_plot), data = data.frame(temp_plot), plotci = FALSE)

      temp1_plot <- RES[[listVARt[k]]]
      temp1_plot <- c(exp(unlist(temp1_plot$estimate)),
                      exp(unlist(temp1_plot$`2.5 %`)),
                      exp(unlist(temp1_plot$`97.5 %`)),
                      exp(unlist(temp1_plot$df)))
      mag_meta_plot$results[1, ] <- temp1_plot
      mag_meta_plot$label <- "All cohorts"

      plT <- .PlotMetaplus2(mag_meta_plot, xlab = "", cex = 1.5, mar = c(2, 5, 4, 2), refline = 1)

      graphics::text(plT[1], 7.8, listVARnames[k], pos = 4, lwd = 2, cex = 2, font = 2)

      temp2_plot <- RES[[j]]$I2

      if (as.numeric(temp2_plot/100) < 0.001) {
        temp3_plot <- bquote(I^2 ~ "<" ~ .(formatC(0.001, digits = 4)))
      }

      if (as.numeric(temp2_plot/100) >= 0.001) {
        temp3_plot <- round(temp2_plot/100, 3)
        temp3_plot <- bquote(I^2 ~ "=" ~ .(formatC(temp3_plot, digits = 4)))
      }

      graphics::text(plT[2], 7.8, temp3_plot, pos = 2, cex = 2)

      j = j + 1
    }
  }
}
