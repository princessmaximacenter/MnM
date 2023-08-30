#' Title
#'
#' @param totalAccDF dataframe with the number of correct and incorrect samples per tumor type size block, e.g. n > 100.
#'
#' @return
#'
#'
#' @examples
getCumulativeAccuracy <- function(totalAccDF) {
  totalAccDF$cumulativeAccuracy  <- 0# Deze maken
  totalAccDF$cumulativeTotal <- 0

  for (i in seq(1:nrow(totalAccDF))) {
    if (i == 1) {
      total <- totalAccDF$total[i]
      totalCorrect <- totalAccDF$totalCorrect[i]
      totalAccDF$cumulativeAccuracy[i] <- totalAccDF$accuracy[i]
      totalAccDF$cumulativeTotal[i] <- total
    } else {
      total <- total + totalAccDF$total[i]
      totalCorrect <- totalCorrect + totalAccDF$totalCorrect[i]
      totalAccDF$cumulativeAccuracy[i] <- round(totalCorrect/total * 100, digits = 1)
      totalAccDF$cumulativeTotal[i] <- total
    }

  }
  return(totalAccDF)
}
