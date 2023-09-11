#' Round all doubles within a vector
#'
#'Function to make sure the combined counts of the normalized RNA-seq data still
#'add up to a million after correction for the ribodepletion protocol.
#'
#' @param x vector with numbers that need to be rounded.
#' @param digits number of decimals the vector needs to be rounded to.
#'
#' @return vector with rounded doubles with _digits_ decimals.
#' @export
#'
smartRound <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}
