#' Boundary of the region covered by two-dimensional data
#'
#' Draws a boundary around the points of two covariates (such as longitude and latitude), following the shape of
#' the data, for the \code{boundary} option of \code{\link{semipred}}: predictions of a surface are shown only
#' inside it, where the data support them.
#'
#' @details The region is the morphological closing of the data with a distance \code{d}. First every place
#' within \code{d} of a data point is taken (a union of discs), which follows any concave shape, then that
#' region is shrunk back by \code{d} less a small \code{buffer}, so that the boundary keeps close to the outer
#' points instead of bulging by \code{d} beyond them. Gaps narrower than about \code{2 d} inside the data are
#' filled; larger empty areas stay outside (holes) and clusters further apart than that become separate pieces.
#' The outline is traced on a grid of \code{n} by \code{n} points with \code{\link[grDevices]{contourLines}}.
#'
#' \code{d} follows the density of the data: by default \code{mult = 2} times the 95th percentile of the distances
#' from each point to its nearest neighbour. A high percentile, not the median, because randomly placed points leave
#' empty patches several times the typical spacing: with twice the median, uniform random designs of 50 to 1000 points
#' gave on average 4 pieces and 12 spurious holes; with the default, one piece and no holes (every point inside), while
#' a real empty area (points on a ring) stays a hole. On the \code{scallop} survey data the default boundary covers
#' 85\% of the same area as the boundary drawn with the original analysis (\code{scallopBdry}; intersection over
#' union) and has every point inside. Larger \code{mult} gives a smoother boundary that bridges larger gaps (large:
#' close to the convex hull); smaller \code{mult} follows the data more tightly and may split it into pieces.
#'
#' Distances are computed on comparable scales: longitude and latitude in degrees, with longitude shrunk by the
#' cosine of the mean latitude (\code{lonlat = TRUE}, the default when the columns are named like longitude
#' and latitude and their values are in range), other covariates after dividing each by its standard deviation.
#'
#' @param x A matrix or data frame with two columns (the two covariates, e.g. longitude and latitude), or a
#'   numeric vector with \code{y}.
#' @param y A numeric vector, when \code{x} is one.
#' @param mult The closing distance \code{d} as a multiple of the 95th percentile of the nearest-neighbour distances
#'   (default 2).
#' @param d The closing distance itself (on the scale described in Details), instead of \code{mult}.
#' @param buffer The margin left around the outer points, as a multiple of the median nearest-neighbour
#'   distance (default 0.5).
#' @param lonlat \code{TRUE} if the columns are longitude and latitude in degrees (in that order), \code{FALSE}
#'   otherwise; \code{NULL} (default) decides from the column names and ranges.
#' @param n The size of the grid the outline is traced on (default 150 by 150).
#' @return A data frame with the two columns of \code{x} (same names), the vertices of the boundary's rings in
#'   order, rings separated by a row of \code{NA} (as \code{\link[graphics]{polygon}} draws them); a point is inside
#'   when it is inside an odd number of rings. Attributes \code{d} (the closing distance) and \code{lonlat}.
#' @seealso \code{\link{semipred}} (\code{boundary = TRUE} uses this function on the model's data).
#' @examples
#' data(scallop); data(scallopBdry)
#' b <- semiboundary(scallop[, c("longitude", "latitude")])
#' plot(scallop$longitude, scallop$latitude, asp = 1.3, pch = 16, cex = 0.6, xlab = "longitude",
#'      ylab = "latitude")
#' polygon(b, border = "blue", lwd = 2)            # from the data
#' polygon(scallopBdry, border = "red")            # drawn with the original analysis
#'
#' \donttest{
#' ## predictions of a surface only where there are data
#' data(scallopKnots)
#' fit <- semireg(asinh(tot.catch) ~ s(longitude, latitude, bs = "ztps",
#'                                    xt = list(knots = as.matrix(scallopKnots))), data = scallop, prt = FALSE)
#' semipred(fit, covariate = c("longitude", "latitude"), boundary = b, prt = FALSE)$plt
#' }
#' @export
semiboundary <- function(x, y = NULL, mult = 2, d = NULL, buffer = 0.5, lonlat = NULL, n = 150) {
  xy <- if (is.null(y)) as.data.frame(x) else data.frame(x = x, y = y)
  if (ncol(xy) != 2 || !all(vapply(xy, is.numeric, logical(1)))) {
    stop("'x' must have two numeric columns (or be a numeric vector with 'y').", call. = FALSE)
  }
  nms <- names(xy)
  xy <- unique(as.matrix(xy[stats::complete.cases(xy), , drop = FALSE]))
  if (nrow(xy) < 3) stop("At least three distinct points are needed for a boundary.", call. = FALSE)
  if (is.null(lonlat)) {
    lonlat <- grepl("^(lon|long|longitude|x_?lon)", nms[1], ignore.case = TRUE) &&
      grepl("^(lat|latitude|y_?lat)", nms[2], ignore.case = TRUE) &&
      all(abs(xy[, 1]) <= 180) && all(abs(xy[, 2]) <= 90)
  }
  sc <- if (lonlat) c(1 / cos(mean(xy[, 2]) * pi / 180), 1) else apply(xy, 2, stats::sd)
  if (any(!is.finite(sc) | sc <= 0)) stop("Each covariate must vary for a boundary.", call. = FALSE)
  z <- sweep(xy, 2, sc, "/")
  nn <- pm_nn_dist(z, z, self = TRUE)
  m <- stats::median(nn)
  if (is.null(d)) d <- mult * stats::quantile(nn, 0.95, names = FALSE)
  if (!is.numeric(d) || length(d) != 1 || !is.finite(d) || d <= 0) stop("'d' must be a positive number.", call. = FALSE)
  shrink <- d - buffer * m
  # grid over the data and a margin of 2 d, so the outline closes inside it
  gx <- seq(min(z[, 1]) - 2 * d, max(z[, 1]) + 2 * d, length.out = n)
  gy <- seq(min(z[, 2]) - 2 * d, max(z[, 2]) + 2 * d, length.out = n)
  G <- as.matrix(expand.grid(gx, gy))
  # dilation: the distance to the data; erosion: the distance from points of
  # the dilated region to its outline
  dist1 <- pm_nn_dist(G, z)
  cl <- grDevices::contourLines(gx, gy, matrix(dist1, n), levels = d)
  rings <- if (shrink <= 0) cl else {
    edge <- do.call(rbind, lapply(cl, function(l) cbind(l$x, l$y)))
    inside <- dist1 <= d
    dist2 <- numeric(nrow(G))
    dist2[inside] <- pm_nn_dist(G[inside, , drop = FALSE], edge)
    grDevices::contourLines(gx, gy, matrix(dist2, n), levels = shrink)
  }
  if (!length(rings)) stop("No boundary could be drawn; try a larger 'mult'.", call. = FALSE)
  out <- do.call(rbind, lapply(seq_along(rings), function(i) {
    r <- cbind(rings[[i]]$x * sc[1], rings[[i]]$y * sc[2])
    if (i < length(rings)) rbind(r, NA) else r
  }))
  out <- stats::setNames(as.data.frame(out), nms)
  structure(out, d = d, lonlat = lonlat)
}

## Distance from each row of A to the nearest row of B (two columns), in
## chunks; self = TRUE: A is B and a point's distance to itself is left out.
pm_nn_dist <- function(A, B, self = FALSE, chunk = 1000) {
  out <- numeric(nrow(A))
  for (s in seq(1, nrow(A), by = chunk)) {
    i <- s:min(nrow(A), s + chunk - 1)
    d2 <- outer(A[i, 1], B[, 1], "-")^2 + outer(A[i, 2], B[, 2], "-")^2
    if (self) d2[cbind(seq_along(i), i)] <- Inf
    out[i] <- sqrt(apply(d2, 1, min))
  }
  out
}
