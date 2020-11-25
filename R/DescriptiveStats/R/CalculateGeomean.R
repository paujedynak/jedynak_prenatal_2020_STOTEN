#' Calculates geometrical mean
#'
#' @param data A dataframe
#'
#' @return A scalar of geometrical mean
#' @export
#' @importFrom Rmisc CI
#' @importFrom stats na.omit
#'
CalculateGeomean <- function(data) {
  geomean <- data %>%
    stats::na.omit() %>%
    log() %>%
    Rmisc::CI(ci = 0.95) %>%
    exp()

  return(geomean)
}
