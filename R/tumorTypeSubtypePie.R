#' Plot pie-chart area with label and color
#'
#' Function to plot part of a piechart, which is required to 'build up' the
#' final piechart layer by layer.
#'
#' @param x Data containing the numbers of a specific tumor (sub)type and their names.
#' @param labels Names for on the pie-chart
#' @param edges ?
#' @param radius What should be the radius of the whole pie-chart?
#' @param clockwise Do we want to go clockwise or counterclock-wise in the pie-chart?
#' @param init.angle At which angle should we start?
#' @param density ?
#' @param angle ?
#' @param col Which colors do you want to use within the pie-chart?
#' @param border Should there be a border? If so, what color should the border have?
#' @param lty Normally what type of line should be applied.
#' @param main What title should the pie-chart get? In our case, we do not use it.
#' @param ... If you want to specify other things.
#'
#' @return Piece of the pie-chart.

tumorTypeSubtypePie <- function(x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE,
                 init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45,
                 col = NULL, border = NULL, lty = NULL, main = NULL, ...)
{
  if (!is.numeric(x) || any(is.na(x) | x < 0))
    stop("'x' values must be positive.")
  if (is.null(labels))
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L])
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  #if (is.null(col))
  #  col <- if (is.null(density))
  #    c("white", "lightblue", "mistyrose", "lightcyan",
  #             "lavender", "cornsilk")
  #             else par("fg")
  if (!is.null(col))
    col <- rep_len(col, nx)
  if (!is.null(border))
    border <- rep_len(border, nx)
  if (!is.null(lty))
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density))
    density <- rep_len(density, nx)
  twopi <- if (clockwise)
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p), an=t2p)
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i],
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      #lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
      text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE,
           srt = ifelse(P$x < 0, P$an/pi*180+180, P$an/pi*180),
           adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}
