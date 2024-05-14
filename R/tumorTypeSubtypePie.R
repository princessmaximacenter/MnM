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
#' @import grDevices

tumorTypeSubtypePie <- function(x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE,
                 init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45,
                 col = NULL, border = NULL, lty = NULL, main = NULL, ...)
{
  if (!base::is.numeric(x) || base::any(base::is.na(x) | x < 0))
    stop("'x' values must be positive.")
  if (base::is.null(labels))
    labels <- base::as.character(base::seq_along(x))
  else labels <- grDevices::as.graphicsAnnot(labels)
  x <- c(0, base::cumsum(x)/base::sum(x))
  dx <- base::diff(x)
  nx <- base::length(dx)
  graphics::plot.new()
  pin <- graphics::par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L])
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  grDevices::dev.hold()
  base::on.exit(grDevices::dev.flush())
  graphics::plot.window(xlim, ylim, "", asp = 1)

  if (!base::is.null(col))
    col <- base::rep_len(col, nx)
  if (!base::is.null(border))
    border <- base::rep_len(border, nx)
  if (!base::is.null(lty))
    lty <- base::rep_len(lty, nx)
  angle <- base::rep(angle, nx)
  if (!base::is.null(density))
    density <- base::rep_len(density, nx)
  twopi <- if (clockwise)
    -2 * base::pi
  else 2 * base::pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * base::pi/180
    list(x = radius * base::cos(t2p), y = radius * base::sin(t2p), an=t2p)
  }
  for (i in 1L:nx) {
    n <- base::max(2, base::floor(edges * dx[i]))
    P <- t2xy(base::seq.int(x[i], x[i + 1], length.out = n))
    graphics::polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i],
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(base::mean(x[i + 0:1]))
    lab <- base::as.character(labels[i])
    if (!base::is.na(lab) && base::nzchar(lab)) {
      graphics::text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE,
           srt = ifelse(P$x < 0, P$an/pi*180+180, P$an/pi*180),
           adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
  graphics::title(main = main, ...)
  base::invisible(NULL)
}
