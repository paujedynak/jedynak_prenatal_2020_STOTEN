#' Function to return an error if there are NAs
#'
#' @param q A name for fit coefficient
#' @param u u
#'
#' @return Error if NAs
#' @importFrom stats coef vcov df.residual
#'
.ExpandVcov <- function(q, u) {
    err <- is.na(q)
    return(u)
    ## if (all(!err)) return(u) k <- length(q) v <- names(q) z <- u for (i in 1:ncol(z)){ if (err[i]) { rbind(z[,],NA,z[,]) j
    ## <- j + 1 up <- } j <- j + 1 z[i,] <- u[j,] z[,i] <- u[,j] }

    ## z <- matrix(NA, ncol=k, nrow=k, dimnames = list(v,v)) idx <- (is.na()) j <- 0 for (i in 1:k){ if (err[i]) next j <- j
    ## + 1 z[i,] <- u[j,] z[,i] <- u[,j] } return(z)
}

#' Execution of rubin's rule depending on a sample size
#'
#' @param m m
#' @param lambda lambda
#' @param dfcom dfcom
#' @param method method
#'
#' @return A number
#'
.MiceDf <- function(m, lambda, dfcom, method) {
    if (is.null(dfcom)) {
        dfcom <- 999999
        warning("Large sample assumed.")
    }
    lambda[lambda < 1e-04] <- 1e-04
    dfold <- (m - 1) / lambda ^ 2
    dfobs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda)
    df <- dfold * dfobs / (dfold + dfobs)
    if (method != "smallsample")
        df <- dfold  ## Rubin 1987, 3.1.6, Van Buuren 2012, 2.30, added 31/10/2012
    return(df)
}

#object=modOne;method = "smallsample"
#' Modified pool function for results of regressions on stacked imputed dataset
#'
#' @param object An object which is a result of a simple regression (one outcome ~ one xposure at a time) run on the mids object created using the modified as.mids function (as.mids.xavi)
#' @param method default = "smallsample"
#'
#' @return A pooled mipo object
#' @import mice
#' @import lme4

.Pool2 <- function(object, method = "smallsample") {

    call <- match.call()

    if (!mice::is.mira(object)) {
        stop("The object must have class 'mira'")
    }

    m <- length(object$analyses)
    fa <- mice::getfit(object, 1)

    if (m == 1) {
        warning("Number of multiple imputations m=1. No pooling done.")
        return(fa)
    }

    analyses <- mice::getfit(object)

    #if (class(fa)[1] == "lme" && !requireNamespace("nlme", quietly = TRUE)) {
    #   stop("Package 'nlme' needed fo this function to work. Please install it.", call. = FALSE)
    #}

    if ((class(fa)[1] == "mer" || class(fa)[1] == "lmerMod" ||
         inherits(fa, "merMod")) && !requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' needed fo this function to work. Please install it.", call. = FALSE)
    }

    mess <- try(stats::coef(fa), silent = TRUE)

    if (inherits(mess, "try-error")) {
        stop("Object has no coef() method.")
    }

    mess <- try(stats::vcov(fa), silent = TRUE)

    if (inherits(mess, "try-error")) {
        stop("Object has no vcov() method.")
    }

    if (class(fa)[1] == "mer" || class(fa)[1] == "lmerMod" ||
        inherits(fa, "merMod")) {
        k <- length(lme4::fixef(fa))
        names <- names(lme4::fixef(fa))
    }

    else if (class(fa)[1] == "polr") {
        k <- length(stats::coef(fa)) + length(fa$zeta)
        names <- c(names(stats::coef(fa)), names(fa$zeta))

    } else {
        fa <- unlist(lapply(object$analyses, function(x) x$coef[!is.na(x$coef)]))
        temp <- table(names(fa))
        names2 <- names(temp)[temp == length(object$analyses)]
        names <- names(object$analyses[[1]]$coef)
        names <- names[names %in% names2]
        k <- length(names)
    }

    qhat <- matrix(NA, nrow = m, ncol = k, dimnames = list(seq_len(m), names))
    u <- array(NA, dim = c(m, k, k), dimnames = list(seq_len(m), names, names))

    for (i in seq_len(m)) {
        fit <- analyses[[i]]

        if (class(fit)[1] == "mer") {
            qhat[i, ] <- lme4::fixef(fit)
            ui <- as.matrix(stats::vcov(fit))

            if (ncol(ui) != ncol(qhat)) {
                stop("Different number of parameters: class mer, fixef(fit): ", ncol(qhat), ", as.matrix(vcov(fit)): ", ncol(ui))
            }

            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }

        else if (class(fit)[1] == "lmerMod" || inherits(fa, "merMod")) {
            qhat[i, ] <- lme4::fixef(fit)
            ui <- stats::vcov(fit)

            if (ncol(ui) != ncol(qhat)) {
                stop("Different number of parameters: class lmerMod, fixed(fit): ", ncol(qhat), ", vcov(fit): ", ncol(ui))
            }

            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }

        else if (class(fit)[1] == "lme") {
            qhat[i, ] <- fit$coefficients$fixed
            ui <- stats::vcov(fit)

            if (ncol(ui) != ncol(qhat)) {
                stop("Different number of parameters: class lme, fit$coefficients$fixef: ", ncol(qhat), ", vcov(fit): ", ncol(ui))
            }

            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }

        else if (class(fit)[1] == "polr") {
            qhat[i, ] <- c(stats::coef(fit), fit$zeta)
            ui <- stats::vcov(fit)

            if (ncol(ui) != ncol(qhat)) {
                stop("Different number of parameters: class polr, c(coef(fit, fit$zeta): ", ncol(qhat), ", vcov(fit): ", ncol(ui))
            }

            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }

        else if (class(fit)[1] == "survreg") {
            qhat[i, ] <- stats::coef(fit)
            ui <- stats::vcov(fit)
            parnames <- dimnames(ui)[[1]]
            select <- !(parnames %in% "Log(scale)")
            ui <- ui[select, select]

            if (ncol(ui) != ncol(qhat)) {
                stop("Different number of parameters: class survreg, coef(fit): ", ncol(qhat), ", vcov(fit): ", ncol(ui))
            }

            u[i, , ] <- array(ui, dim = c(1, dim(ui)))

        } else {
            qhat[i, ] <- stats::coef(fit)[names(stats::coef(fit)) %in% names]
            ui <- stats::vcov(fit)[colnames(stats::vcov(fit)) %in% names, rownames(stats::vcov(fit)) %in% names]
            ui <- .ExpandVcov(qhat[i, ], ui)
            #            if (ncol(ui) != ncol(qhat))
            #                stop("Different number of parameters: coef(fit): ",
            #                  ncol(qhat), ", vcov(fit): ", ncol(ui))
            u[i, , ] <- array(ui, dim = c(1, dim(ui)))
        }
    }

    qbar <- apply(qhat, 2, mean)
    ubar <- apply(u, c(2, 3), mean)
    e <- qhat - matrix(qbar, nrow = m, ncol = k, byrow = TRUE)
    b <- (t(e) %*% e) / (m - 1)
    t <- ubar + (1 + 1 / m) * b
    r <- (1 + 1 / m) * diag(b/ubar)
    lambda <- (1 + 1 / m) * diag(b / t)
    dfcom <- stats::df.residual(object)
    df <- .MiceDf(m, lambda, dfcom, method)
    fmi <- (r + 2 / (df + 3)) / (r + 1)
    names(r) <- names(df) <- names(fmi) <- names(lambda) <- names
    fit <- list(call = call, call1 = object$call, call2 = object$call1,
                nmis = object$nmis, m = m, qhat = qhat, u = u, qbar = qbar,
                ubar = ubar, b = b, t = t, r = r, dfcom = dfcom, df = df,
                fmi = fmi, lambda = lambda)
    oldClass(fit) <- c("mipo", oldClass(object))
    return(fit)
}


#' Function
#'
#' @param est est
#' @param se se
#' @param level level
#'
#' @return A dataframe with estimates and CIs for each stratum
#'
.Modif <- function(est = est, se = se, level) {
    # wh <- as.numeric(as.character(gsub("Cohort","",sapply(strsplit(names(est),"[:]"),function(x) x[1]))))
    # if(!all(names(est)==colnames(se))) stop()
    # print(sqrt(abs(se)))
    wh <- grepl("[:]", names(est))
    est[wh] <- est[wh] + est[!wh]
    se2 <- diag(se)
    se2[wh] <- se2[wh] + se2[!wh] + 2 * se[wh, !wh]
    res <- cbind(est = est,se = sqrt(se2))
    rownames(res) <- level
    data.frame(res)
}




