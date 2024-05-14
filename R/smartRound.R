#' Round all doubles within a vector
#'
#'Function to make sure the combined counts of the normalized RNA-seq data still
#'add up to a million after correction for the ribodepletion protocol.
#'
#' @param x vector with numbers that need to be rounded.
#' @param digits number of decimals the vector needs to be rounded to.
#'
#' @return vector with rounded doubles with _digits_ decimals.
#' @import utils
#'
smartRound <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- base::floor(x)
  indices <- utils::tail(order(x-y), base::round(base::sum(x)) - base::sum(y))
  y[indices] <- y[indices] + 1
  y / up
}
