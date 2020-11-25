#' Correction of the as.mids function integrated in the mice package
#'
#' @param data A complete stacked dataset (mice object) containing outcome and exposures of interest, covariate variables, .id, .imp,
#' @param .imp Position of the .imp column
#' @param .id Position of the .id column
#'
#' @return A modified mids object
#' @import mice

.AsMidsXavi <- function(data, .imp, .id) {

  starting_values <- function(dat, m) {
    ini <- mice::mice(dat, m = m, maxit = 0)
    dat_comp <- mice::complete(ini, 1)

    for (i.col in 1:dim(dat)[2]) {

      if (sum(is.na(dat_comp[, i.col])) > 0) {
        # variables with missing values in the original dataset, without initial
        # values (imputations) in the completed dataset
        y <- dat_comp[, i.col]
        ry <- !is.na(y)
        mat_aux <- matrix(NA,nrow = sum(is.na(y)), ncol = m)

        for (i.imp in seq_len(m)) {
          mat_aux[, i.imp] <- mice::mice.impute.sample(y, ry)
        }

        df <- data.frame(mat_aux)
        rownames(df) <- rownames(dat_comp)[is.na(y)]
        names(df) <- as.character(seq_len(m))

        ini$imp[[i.col]] <- df
      }
    }

    return(ini)
  }

  .AsMidsMod <- function(data, .imp = 1, .id = 2) {

    data[, .imp] <- as.numeric(as.character(data[, .imp]))
    data[, .id] <- as.numeric(as.character(data[, .id]))
    data <- data[order(data[, .imp]),]
    imps <- as.numeric(as.character(data[, .imp]))

    if (is.factor(imps)) {
      m <- max(as.numeric(levels(imps))[imps])

    } else {
      m <- max(imps)
    }

    if (all(imps != 0)) {
      stop("the non-imputed values need to be provided in data")
    }

    ini <- starting_values(data[imps == 0, ], m = m)
    names <- names(ini$imp)

    if (!is.null(.id)) {
      rownames(ini$data) <- data[imps == 0, .id]
    }

    for (i in seq_along(names)) {

      for (j in seq_len(m)) {

        if (!is.null(ini$imp[[names[i]]])) {

          indic <- imps == j & is.na(data[imps == 0, names[i]])
          ini$imp[[names[i]]][j] <- data[indic, names[i]]
        }
      }
    }

    return(ini)
  }

  result <- .AsMidsMod(data, .imp, .id)

  return(result)
}
