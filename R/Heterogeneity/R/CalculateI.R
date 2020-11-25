# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c(".imp", ".id", "var"))

#' # Calculate heterogeneity measure (I squared) - CODE and COMMENTS come from Xavier Basagaña
#'
#' @param data A complete stacked dataset (mice object), including .imp = 0
#' @param outcome A string defining an outcome ("h_sdq_external" or "h_sdq_ternal")
#' @param conf_list A vector containing all confounder variable names
#' @param selected_exposures A vector containing all exposure variable names
#' @param strat_variable A string describing stratifying variable (normally "cohort")
#'
#' @return A list containing RES (estimates and CIs for all strata) and RESv (estimates and CIs for each strata separately)
#' @import mice
#' @import dplyr
#' @importFrom stats confint
#' @importFrom MASS glm.nb
#' @importFrom MAc mareg

.CalculateI <- function(data, outcome, conf_list, selected_exposures, strat_variable) {

  ### formula and restricted data
  form <- paste(outcome, "~ 1")

  if (length(conf_list) > 0) {

    for (item in conf_list) {

      form <- paste(form, "+", item)
    }
  }

  dataT <- dplyr::select(data, outcome, selected_exposures, .imp, .id, conf_list, strat_variable)

  ### linear models testing exposure, and exposure- strat_variable interaction, with I2
  imp.completed <- .AsMidsXavi(data = dataT, .imp = which(colnames(dataT) == ".imp"), .id = which(colnames(dataT) == ".id"))

  RESv <- RES <- list()

  for (i in seq_along(selected_exposures)) {
    print(paste0(i, "/", length(selected_exposures)))

    ## model with exposure
    temp <- summary(mice::pool(with(imp.completed, MASS::glm.nb(as.formula(paste(form, "+", selected_exposures[i]))))), conf.int = TRUE)

    ## model with exposure-strat_variable interaction
    RES[[i]] <- temp[nrow(temp),]

    modOneAll <- .Pool2(with(imp.completed, MASS::glm.nb(as.formula(paste(form, "+", strat_variable, "*", selected_exposures[i])))))

    wh <- c(which(rownames(modOneAll$t) == selected_exposures[i]), grep(paste(":", selected_exposures[i], sep = ""), rownames(modOneAll$t)))

    # The following section assumes the following order of the cohorts: "BiB"  "EDEN" "INMA" "KANC" "RHEA" and this order must be kept!
    wh2 <- c("BiB", gsub(strat_variable, "", sapply(strsplit(names(modOneAll$qbar)[wh[-1]], "[:]"), function(x) x[1])))

    temp <-  .Modif(modOneAll$qbar[wh], modOneAll$t[wh, wh], level = wh2)
    RESv[[i]] <- temp
    temp$var <- temp$se ^ 2
    m0 <- MAc::mareg(est ~ 1, var = var, data = data.frame(temp))
    RES[[i]] <- c(RES[[i]], I2 = stats::confint(m0)$random[3, 1])
  }

  names(RESv) <- names(RES) <- selected_exposures
  obj <- list("RES" = RES, "RESv" = RESv)

  return(obj)
}
