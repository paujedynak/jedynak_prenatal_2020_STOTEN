#' Plot heterogeneity (forest plots)
#'
#' @param x A dataset of mag_meta format (for all strata, e.g. All cohorts)
#' @param ... Parameters of the metafor::forest plot (cex, margins etc.)
#' @param extrameta default NULL
#' @param mar Default margins
#'
#' @return A forest plot
#' @import metafor
#' @import graphics

.PlotMetaplus2 <- function(x, ..., extrameta = NULL, mar = c(1, 1, 1, 1)) {
  #x=mag.meta;xlab = listVARt[i];  cex = 0.75;mar=c(0.5,4,0,0.5)

  if (!inherits(x, "metaplus")) {
    stop("Use only with 'metaplus' xs.\n")
  }

  if (!is.null(extrameta)) {
    if (any(!sapply(extrameta, inherits, what = "metaplus"))) {
      stop("Use only with 'metaplus' xs.\n")
    }
  }

  if (x$justfit) {
    stop("Cannot use with objects fitted with justfit = TRUE")
  }

  if (is.null(extrameta)) {
    nsummaries <- 1

    } else {
      nsummaries <- length(extrameta) + 1
    }

  graphics::par(mar = mar)

  temp <- metafor::forest(x = x$yi, sei = x$sei, slab = x$slab, ylim = c(-0.5 - nsummaries, length(x$yi) + 3),
                          mar = c(0.5, 4, 0, 0.5), ...)
  sumxi <- x$results[1, 1]
  sumci_lb <- x$results[1, 2]
  sumci_ub <- x$results[1, 3]
  mlab <- x$label

  if (!is.null(extrameta)) {
    for (iadd in seq_along(extrameta)) {
      thex <- extrameta[[iadd]]
      sumxi <- c(sumxi, thex$results[1, 1])
      sumci_lb <- c(sumci_lb,thex$results[1, 2])
      sumci_ub <- c(sumci_ub,thex$results[1, 3])
      mlab <- c(mlab, thex$label)
    }
  }

  extravars <- list(...)

  if (length(extravars) == 0) {
    metafor::addpoly.rma(x = sumxi, ci.lb = sumci_lb, ci.ub = sumci_ub, mlab = mlab)

  } else {
    # only pass relevent parameters to addpoly
    extravars <- list(transf = extravars$transf, atransf = extravars$atransf, targs = extravars$targs,
                      efac = extravars$efac, cex = extravars$cex, digits = extravars$digits)
    extravars <- extravars[!sapply(extravars, is.null)]
    extravars <- c(list(x = sumxi, ci.lb = sumci_lb, ci.ub = sumci_ub, mlab = mlab), extravars)
    do.call("addpoly", extravars)
  }

  graphics::abline(h = 0)

  return(temp$xlim)
}
