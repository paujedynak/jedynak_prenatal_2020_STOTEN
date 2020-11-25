#' pool_compare2 (some function of Xavi that I do not use) calculating heterogeneity for not small sample size
#'
#' @param fit1 fit1
#' @param fit0 fit0
#' @param data data
#' @param method method
#'
#' @return statistic
#' @importFrom stats model.matrix model.frame formula pf

.PoolCompare2 <- function(fit1, fit0, data = NULL, method = "Wald") {

  .LLlogistic <- function(formula, data, coefs) {

    .logistic <- function(mu) {
      res <- exp(mu) / (1 + exp(mu))

      return(res)
    }

    Xb <- stats::model.matrix(formula, data) %*% coefs
    y <- stats::model.frame(formula, data)[1][, 1]
    p <- .logistic(Xb)
    y <- (y - min(y)) / (max(y) - min(y))
    term1 <- term2 <- rep(0, length(y))
    term1[y != 0] <- y[y != 0] * log(y[y != 0]/p[y != 0])
    term2[y == 0] <- (1 - y[y == 0]) * log((1 - y[y == 0]) / (1 - p[y == 0]))

    return(-(2 * sum(term1 + term2)))
  }

  call <- match.call()
  meth <- match.arg(tolower(method), c("wald", "likelihood"))

  if (!is.mira(fit1) || !is.mira(fit0)) {
    stop("fit1 and fit0 should both have class 'mira'.\n")
  }

  m1 <- length(fit1$analyses)
  m0 <- length(fit0$analyses)

  if (m1 != m0) {
    stop("Number of imputations differs between fit1 and fit0.\n")
  }

  if (m1 < 2) {
    stop("At least two imputations are needed for pooling.\n")
  }

  m <- m1
  est1 <- .Pool2(fit1)
  est0 <- .Pool2(fit0)
  dimQ1 <- length(est1$qbar)
  dimQ2 <- dimQ1 - length(est0$qbar)
  formula1 <- stats::formula(fit1$analyses[[1]])
  formula0 <- stats::formula(fit0$analyses[[1]])
  vars1 <- names(est1$qbar)
  vars0 <- names(est0$qbar)

  if (is.null(vars1) || is.null(vars0)) {
    stop("coefficients do not have names")
  }

  if (dimQ2 < 1) {
    stop("The larger model should be specified first and must be strictly larger than the smaller model.\n")
  }

  if (!setequal(vars0, intersect(vars0, vars1))) {
    stop("The smaller model should be fully contained in the larger model. \n")
  }

  if (meth == "wald") {
    Q <- diag(rep.int(1, dimQ1), ncol = dimQ1)
    where_new_vars = which(!(vars1 %in% vars0))
    Q <- Q[where_new_vars, , drop = FALSE]
    qbar <- Q %*% est1$qbar
    Ubar <- Q %*% est1$ubar %*% (t(Q))
    Bm <- Q %*% est1$b %*% (t(Q))
    rm <- (1 + 1 / m) * sum(diag(Bm %*% (solve(Ubar)))) / dimQ2
    Dm <- (t(qbar)) %*% (solve(Ubar)) %*% qbar/(dimQ2 * (1 + rm))
  }

  if (meth == "likelihood") {
    if (is.null(data)) {
      stop("For method = likelihood the imputed data set (a mids object) should be included.\n")
    }

    devM <- devL <- 0

    for (i in seq_len(m)) {
      devL <- devL + .LLlogistic(formula1, complete(data, i), est1$qbar) -
        .LLlogistic(formula0, complete(data, i), est0$qbar)

      devM <- devM + .LLlogistic(formula1, complete(data, i), est1$qhat[i, ]) -
        .LLlogistic(formula0, complete(data, i), est0$qhat[i, ])
    }

    devL <- devL / m
    devM <- devM / m
    rm <- ((m + 1) / (dimQ2 * (m - 1))) * (devM - devL)
    Dm <- devL/(dimQ2 * (1 + rm))
  }

  v <- dimQ2 * (m - 1)

  if (v > 4) {
    w <- 4 + (v - 4) * ((1 + (1 - 2/v) * (1/rm))^2)

  } else {
    w <- v * (1 + 1/dimQ2) * ((1 + 1/rm)^2)/2

  }

  statistic <- list(call = call, call11 = fit1$call, call12 = fit1$call1,
                    call01 = fit0$call, call02 = fit0$call1, method = method,
                    nmis = fit1$nmis, m = m, qhat1 = est1$qhat, qhat0 = est0$qhat,
                    qbar1 = est1$qbar, qbar0 = est0$qbar, ubar1 = est1$ubar,
                    ubar0 = est0$ubar, Dm = Dm, rm = rm, df1 = dimQ2, df2 = w,
                    pvalue = 1 - stats::pf(Dm, dimQ2, w))

  return(statistic)
}
