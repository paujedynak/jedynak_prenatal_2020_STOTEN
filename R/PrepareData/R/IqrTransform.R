# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c(".imp", "HelixID", "h_sdq_internal", "smoking", "sex", "creatinine", "fat", "cohort"))

#' IQR standardize exposure variables concentrations
#'
#' @param dataset An imputed dataset (mids object or a dataframe)
#'
#' @return A dataframe with IQR transformed variables replacing original ones
#' @export
#' @import mice
#' @import dplyr
#' @importFrom dplyr %>%
#' @importFrom stats quantile

IqrTransform <- function(dataset) {

    if (!(class(dataset) == "mids" | class(dataset) == "data.frame")) {
        stop("data must be a mids object or a data.frame!")
    }

    # Create a function to iqr standardize variables
    calculate_iqr <- function(dataset) {

        # For imputed 'mids' dataset
        if (class(dataset) == "mids") {

            # Create complete stacked imputed dataset, including original dataset
            imp_long <- mice::complete(dataset, action = "long", include = TRUE)

            # choose variables to be transformed (do not transform covariates, sdqs, .imp and .id)
            # Create columns that will contain IQR standardized values
            vars_to_iqr <- iqr <- imp_long %>%
                dplyr::select(.imp, HelixID, ends_with("_log2"), ends_with("_ln"))

            # For each variable to be iqr transformed (except .imp)
            for (i in 3:(ncol(vars_to_iqr))) {

                # calculate iqr for the first imputed dataset
                iqr_first <- stats::quantile(vars_to_iqr[vars_to_iqr$.imp ==
                  1, i], 0.75) - quantile(vars_to_iqr[vars_to_iqr$.imp ==
                  1, i], 0.25)

                # Divide each untransformed value by the iqr value obtained in the previous step
                iqr[, i] <- iqr[, i]/iqr_first
            }

            # Change names of the new iqr transformed variables
            colnames(iqr)[-c(1:2)] <- paste(colnames(iqr)[-c(1:2)],
                "iqr", sep = "_")

            # Merge iqr transformed variables with the rest of the dataset
            transf_data <- merge(dplyr::select(imp_long, .imp:h_sdq_internal, smoking:sex, creatinine, fat), iqr, by = c(".imp", "HelixID"))

        # For non-imputed datasets, other than 'mids'
        } else {

            #choose variables to be transformed (do not transform covariates, sdqs, .imp and .id)
            vars_to_iqr <- iqr <- dplyr::select(dataset, HelixID, ends_with("_log2"),
                ends_with("_ln"))

            for (i in 2:ncol(vars_to_iqr)) {

                # iqr of first imputed dataset
                iqr_first <- quantile(vars_to_iqr[, i], 0.75,
                  na.rm = TRUE) - quantile(vars_to_iqr[, i],
                  0.25, na.rm = TRUE)

                iqr[, i] <- iqr[, i]/iqr_first
            }

            # Change names of the new iqr transformed variables
            colnames(iqr)[-1] <- paste(colnames(iqr)[-1], "iqr", sep = "_")

            # Merge iqr transformed variables with the rest of the dataset
            transf_data <- merge(dplyr::select(dataset, HelixID, cohort:h_sdq_internal, smoking:sex), iqr, by = "HelixID")
        }

        return(transf_data)
    }

    # Obtain iqr standardized values
    iqr_transformed <- calculate_iqr(dataset)

    # For imputed 'mids' dataset
    if (class(dataset) == "mids") {

        # Drop unnecessary levels before closing in the mids object
        iqr_transformed <- droplevels(iqr_transformed) %>%

            # Change back the long stacked dataset into a 'mids' object
            mice::as.mids()
    }

    return(iqr_transformed)
}
