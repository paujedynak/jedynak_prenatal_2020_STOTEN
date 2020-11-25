#' Adjust exposure variables concentrations with creatinine concentration or total fat percentage
#'
#' @param data A dataframe containing variables to be adjusted
#' @param list_fat_adj A character vector with names of variables to be adjusted for total fat percentage
#' @param list_creat_adj A character vector with names of variables to be adjusted for creatinine
#'
#' @return A dataframe containing original and adjusted variables
#' @import dplyr
#' @importFrom tidyselect all_of

.AdjustFatCreat <- function(data, list_fat_adj, list_creat_adj) {

    if (class(data) != "data.frame") {
        stop("data must be a data.frame!")
    }

    if (class(list_fat_adj) != "character" | class(list_creat_adj) != "character") {
        stop("list_fat_adj and list_creat_adj must be a character vector")
    }

    # Adjust for fat
    fat_adj <- sweep((dplyr::select(data, tidyselect::all_of(list_fat_adj)) / 2), 1,
                     (data$fat * 10), "/")

    # Change names of new adjusted variables
    names(fat_adj) <- gsub(x = names(fat_adj), pattern = "_m\\>",
                           replacement = "_adj")

    # Adjust for creatinine
    creat_adj <- sweep(dplyr::select(data, tidyselect::all_of(list_creat_adj)), 1,
                       data$creatinine, "/")

    # Change names of new adjusted variables
    names(creat_adj) <- gsub(x = names(creat_adj), pattern = "_m\\>",
                             replacement = "_adj")

    # Merge adjusted variables with the rest of the dataset
    adj_data <- cbind(data, fat_adj, creat_adj)

    return(adj_data)
}
