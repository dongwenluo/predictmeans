#' Boundary of the region covered by two-dimensional data
#'
#' Draws a boundary around the points of two covariates (such as longitude and latitude, or the coordinates of
#' drill holes), following the shape of the data, for the \code{boundary} option of \code{\link{semipred}}:
#' predictions of a surface are then shown only where the data support them, not in the empty corners of the
#' rectangle the grid covers.
#'
#' @details The region is the morphological closing of the data with a distance \code{d}. First every place
#' within \code{d} of a data point is taken (a union of discs), which follows any concave shape; then that
#' region is shrunk back by \code{d} less a small \code{buffer}, so that the boundary keeps close to the outer
#' points instead of bulging by \code{d} beyond them. Gaps narrower than about \code{2 d} inside the data are
#' filled; larger empty areas stay outside (holes), and clusters further apart than that become separate pieces.
#' The outline is traced on a grid of \code{n} by \code{n} points with \code{\link[grDevices]{contourLines}}
#' (base R: no other package is needed).
#'
#' \code{d} follows the density of the data: by default \code{mult = 2} times the 95th percentile of the distances
#' from each point to its nearest neighbour. A high percentile, not the median, because randomly placed points leave
#' empty patches several times the typical spacing: with twice the median, uniform random designs of 50 to 1000
#' points gave on average 4 pieces and 12 spurious holes; with the default, one piece and no holes (every point
#' inside), while a real empty area (points on a ring) stays a hole. On the \code{scallop} survey the default
#' boundary covers 85\% of the same area as the boundary drawn with the original analysis (\code{scallopBdry};
#' intersection over union), with every point inside.
#'
#' Distances need the two covariates on comparable scales (\code{scale}): \code{"lonlat"}, longitude and latitude
#' in degrees, longitude shrunk by the cosine of the mean latitude (a degree east is shorter than a degree north
#' away from the equator); \code{"same"}, coordinates in the same unit (metres, a local grid), used as they are;
#' \code{"sd"}, two different variables, each divided by its standard deviation. \code{"auto"} (default) takes
#' \code{"lonlat"} when the columns are named like longitude and latitude (in that order) with values in range,
#' else \code{"sd"}.
#'
#' @section Suggestions:
#' \itemize{
#'   \item Look at the boundary with the data before using it (\code{plot()} of the points and \code{polygon()} of
#'     the boundary, as in the examples). It is a choice about where predictions are believable, not an estimate.
#'   \item \code{mult}: the default suits most designs, regular grids and random sampling alike. Use 3 or 4 for
#'     sparse or uneven sampling, or to keep one region where the data come in clusters; use 1 to 1.5 to follow
#'     the corners of a regular grid tightly. Separate pieces for clusters far apart are honest: there are no data
#'     between them.
#'   \item \code{scale}: for map coordinates in metres (UTM, a local grid) give \code{scale = "same"}, so that the
#'     shape is not stretched; for two different covariates (for example temperature and pressure) the default
#'     \code{"sd"} keeps predictions to the combinations that were observed, which avoids extrapolating a surface
#'     into corners without data.
#'   \item The boundary only chooses which grid points \code{semipred()} shows; it does not change the fit. Near
#'     the edge the bands are wider, as they should be.
#'   \item A boundary of your own (a coastline, a farm or field outline) can be given to \code{semipred()} directly:
#'     two columns, several rings separated by a row of \code{NA}.
#'   \item Holes: \code{semipred()} leaves them out (a point is inside when it is inside an odd number of rings).
#'     To shade a region with holes use \code{\link[graphics]{polypath}(b, rule = "evenodd")}; to draw the outline
#'     with ggplot2 number the rings, \code{b$ring <- cumsum(is.na(b[[1]]))}, and use \code{geom_path(aes(group =
#'     ring))}, as in the examples.
#'   \item Large data: about 0.1 s for 150 points, 4 s for 10,000 and 9 s for 50,000; the grid the outline is
#'     traced on gets finer as the data get denser (up to a million cells).
#' }
#'
#' @param x A matrix or data frame with two columns (the two covariates, e.g. longitude and latitude), or a
#'   numeric vector with \code{y}.
#' @param y A numeric vector, when \code{x} is one.
#' @param mult The closing distance \code{d} as a multiple of the 95th percentile of the nearest-neighbour distances
#'   (default 2).
#' @param d The closing distance itself (on the scale of \code{scale}), instead of \code{mult}.
#' @param buffer The margin left around the outer points, as a multiple of the median nearest-neighbour
#'   distance (default 0.5).
#' @param scale How the two covariates are put on comparable scales: \code{"auto"} (default), \code{"lonlat"},
#'   \code{"same"} or \code{"sd"} (see Details).
#' @param min_hole Holes smaller than the area of this many data points (the region's area divided by the number of
#'   points, times \code{min_hole}; default 20) are filled: the empty patches that random sampling leaves, not
#'   areas without data. \code{0} keeps every hole.
#' @param n The least number of grid cells along the longer side of the data the outline is traced on (default
#'   150); finer where the data are dense (cells of at most \code{d / 4}, up to a million cells).
#' @return A data frame with the two columns of \code{x} (same names), the vertices of the boundary's rings in
#'   order, rings separated by a row of \code{NA} (as \code{\link[graphics]{polygon}} draws them); a point is inside
#'   when it is inside an odd number of rings. Attributes \code{d} (the closing distance) and \code{scale}.
#' @seealso \code{\link{semipred}} (\code{boundary = TRUE} uses this function on the model's data), the data sets
#'   \code{\link{scallop}} and \code{\link{copper}}.
#' @examples
#' library(predictmeans)
#' ## 1. The scallop survey (longitude and latitude): the boundary from the data (blue) and the
#' ##    boundary drawn with the original analysis (red)
#' data(scallop); data(scallopBdry)
#' b <- semiboundary(scallop[, c("longitude", "latitude")])
#' attr(b, "scale")                                   # "lonlat", from the column names
#' plot(latitude ~ longitude, data = scallop, asp = 1/cos(40*pi/180), pch = 16, cex = 0.6,
#'      xlim = range(b$longitude, na.rm = TRUE), ylim = range(b$latitude, na.rm = TRUE))
#' polygon(b, border = "blue", lwd = 2)
#' polygon(scallopBdry, border = "red", lwd = 2)
#' legend("bottomright", c("semiboundary()", "scallopBdry"), col = c("blue", "red"), lwd = 2, bty = "n")
#'
#' ## 2. Copper drill holes (coordinates in the same unit): how 'mult' changes the region;
#' ##    with mult = 1 the hole far to the left is a separate piece
#' data(copper)
#' holes <- unique(copper[, c("xcoord", "ycoord")])
#' b1 <- semiboundary(holes, mult = 1, scale = "same")
#' b2 <- semiboundary(holes, scale = "same")
#' b4 <- semiboundary(holes, mult = 4, scale = "same")
#' plot(holes, asp = 1, pch = 16, xlim = range(b4$xcoord, na.rm = TRUE),
#'      ylim = range(b4$ycoord, na.rm = TRUE))
#' polygon(b1, border = "orange", lwd = 2)
#' polygon(b2, border = "blue", lwd = 2)
#' polygon(b4, border = "darkgreen", lwd = 2)
#' legend("topleft", c("mult = 1", "mult = 2 (default)", "mult = 4"),
#'        col = c("orange", "blue", "darkgreen"), lwd = 2, bty = "n")
#'
#' ## 3. Pieces and holes: random points on a ring (the empty middle is a hole) and a separate
#' ##    cluster; polypath() shades the region with the hole left open
#' set.seed(1)
#' th <- runif(400, 0, 2*pi); r <- runif(400, 2, 3)
#' pts <- data.frame(x = c(r*cos(th), runif(80, -3, -2)), y = c(r*sin(th), runif(80, 3.5, 4.5)))
#' br <- semiboundary(pts)
#' plot(pts, asp = 1, pch = 16, cex = 0.4)
#' polypath(br, col = adjustcolor("skyblue", 0.4), border = "blue", lwd = 2, rule = "evenodd")
#' sum(is.na(br$x)) + 1                               # 3 rings: the ring's two edges and the cluster
#'
#' \donttest{
#' ## Surfaces shown only inside the boundary, with its outline and the data
#' library(ggplot2)
#' outline <- function(b) { b$ring <- cumsum(is.na(b[[1]])); b }
#' data(scallopKnots)
#' fs <- semireg(asinh(tot.catch) ~ s(longitude, latitude, bs = "ztps",
#'                                   xt = list(knots = as.matrix(scallopKnots))), data = scallop, prt = FALSE)
#' semipred(fs, covariate = c("longitude", "latitude"), boundary = b, trans = sinh)$plt +
#'   geom_path(data = outline(b), aes(longitude, latitude, group = ring), inherit.aes = FALSE,
#'             na.rm = TRUE)
#'
#' data(copperKnots)
#' fc <- semireg(grade ~ s(xcoord, ycoord, bs = "ztps", xt = list(knots = as.matrix(copperKnots))) +
#'                 s(zcoord, bs = "os"), data = copper, prt = FALSE)
#' semipred(fc, covariate = c("xcoord", "ycoord"), boundary = b2)$plt +
#'   geom_path(data = outline(b2), aes(xcoord, ycoord, group = ring), inherit.aes = FALSE, na.rm = TRUE) +
#'   geom_point(data = holes, aes(xcoord, ycoord), inherit.aes = FALSE, size = 0.8)
#' # boundary = TRUE draws the boundary from the model's data (scale "auto": here "sd")
#' semipred(fc, covariate = c("xcoord", "ycoord"), boundary = TRUE)$plt
#' }
#' @export
semiboundary <- function(x, y = NULL, mult = 2, d = NULL, buffer = 0.5, scale = c("auto", "lonlat", "same", "sd"),
                         min_hole = 20, n = 150) {
  scale <- match.arg(scale)
  xy <- if (is.null(y)) as.data.frame(x) else data.frame(x = x, y = y)
  if (ncol(xy) != 2 || !all(vapply(xy, is.numeric, logical(1)))) {
    stop("'x' must have two numeric columns (or be a numeric vector with 'y').", call. = FALSE)
  }
  nms <- names(xy)
  xy <- unique(as.matrix(xy[stats::complete.cases(xy), , drop = FALSE]))
  if (nrow(xy) < 3) stop("At least three distinct points are needed for a boundary.", call. = FALSE)
  sds <- apply(xy, 2, stats::sd)
  if (any(!is.finite(sds) | sds <= 0)) stop("Each covariate must vary for a boundary.", call. = FALSE)
  if (scale == "auto") {
    ll <- grepl("^(lon|long|longitude|x_?lon)", nms[1], ignore.case = TRUE) &&
      grepl("^(lat|latitude|y_?lat)", nms[2], ignore.case = TRUE) &&
      all(abs(xy[, 1]) <= 180) && all(abs(xy[, 2]) <= 90)
    scale <- if (ll) "lonlat" else "sd"
  }
  sc <- switch(scale, lonlat = c(1 / cos(mean(xy[, 2]) * pi / 180), 1), same = c(1, 1), sd = sds)
  z <- sweep(xy, 2, sc, "/")
  # nearest-neighbour distances (of up to 2000 points spread through the data,
  # each against all points: their quantiles, at a fraction of the cost)
  q <- if (nrow(z) <= 2000) seq_len(nrow(z)) else unique(round(seq(1, nrow(z), length.out = 2000)))
  nn <- pm_nn_dist(z[q, , drop = FALSE], z, self = q)
  m <- stats::median(nn)
  if (is.null(d)) d <- mult * stats::quantile(nn, 0.95, names = FALSE)
  if (!is.numeric(d) || length(d) != 1 || !is.finite(d) || d <= 0) stop("'d' must be a positive number.", call. = FALSE)
  # grid over the data and a margin of 2 d (so the outline closes inside it),
  # the same spacing on both axes: at least n cells along the longer side and
  # at most d / 4 (the outline then lies within about 1/8 of d), at most
  # 10^6 cells
  lo <- apply(z, 2, min) - 2 * d
  span <- apply(z, 2, max) + 2 * d - lo
  h <- min(max(span) / (n - 1), d / 4)
  h <- max(h, sqrt(prod(span) / 1e6))
  shrink <- d - max(buffer * m, h)                 # a margin of at least one cell
  gx <- seq(lo[1], lo[1] + span[1] + h, by = h)
  gy <- seq(lo[2], lo[2] + span[2] + h, by = h)
  # dilation: the distance to the data; erosion: the distance from places in
  # the dilated region to its outline (both exact up to the distances that
  # matter, pm_grid_dist())
  dist1 <- pm_grid_dist(gx, gy, z, cap = 1.5 * d)
  cl <- grDevices::contourLines(gx, gy, dist1, levels = d)
  rings <- if (shrink <= 0) cl else {
    edge <- do.call(rbind, lapply(cl, function(l) cbind(l$x, l$y)))
    dist2 <- pm_grid_dist(gx, gy, edge, cap = d)
    dist2[dist1 > d] <- 0
    grDevices::contourLines(gx, gy, dist2, levels = shrink)
  }
  if (!length(rings)) stop("No boundary could be drawn; try a larger 'mult'.", call. = FALSE)
  # holes (rings inside an odd number of others) smaller than the area of
  # min_hole data points are filled: empty patches of random sampling, not
  # areas without data
  if (length(rings) > 1) {
    rm <- lapply(rings, function(l) cbind(l$x, l$y))
    area <- vapply(rm, function(r) abs(sum(r[, 1] * c(r[-1, 2], r[1, 2]) - c(r[-1, 1], r[1, 1]) * r[, 2])) / 2,
                   numeric(1))
    hole <- vapply(seq_along(rm), function(i) {
      sum(vapply(rm[-i], function(o) pm_points_in_poly(rm[[i]][1, , drop = FALSE], o), logical(1))) %% 2 == 1
    }, logical(1))
    per_point <- (sum(area[!hole]) - sum(area[hole])) / nrow(z)
    rings <- rings[!(hole & area < min_hole * per_point)]
  }
  out <- do.call(rbind, lapply(seq_along(rings), function(i) {
    r <- cbind(rings[[i]]$x * sc[1], rings[[i]]$y * sc[2])
    if (i < length(rings)) rbind(r, NA) else r
  }))
  out <- stats::setNames(as.data.frame(out), nms)
  structure(out, d = d, scale = scale)
}

## Distance from each row of A to the nearest row of B (two columns), in
## chunks; self: the rows of B that the rows of A are (left out).
pm_nn_dist <- function(A, B, self = NULL) {
  out <- numeric(nrow(A))
  chunk <- max(1L, floor(4e6 / nrow(B)))         # about 32 MB of distances at a time
  for (s in seq(1, nrow(A), by = chunk)) {
    i <- s:min(nrow(A), s + chunk - 1)
    d2 <- outer(A[i, 1], B[, 1], "-")^2 + outer(A[i, 2], B[, 2], "-")^2
    if (!is.null(self)) d2[cbind(seq_along(i), self[i])] <- Inf
    out[i] <- sqrt(d2[cbind(seq_along(i), max.col(-d2, ties.method = "first"))])
  }
  out
}

## Distance from each point of the grid gx x gy to the nearest point of z,
## exact up to 'cap' (larger distances are returned as cap): each point
## updates only the grid cells within cap of it.
pm_grid_dist <- function(gx, gy, z, cap) {
  nx <- length(gx)
  ny <- length(gy)
  dx <- gx[2] - gx[1]
  dy <- gy[2] - gy[1]
  wx <- ceiling(cap / dx)
  wy <- ceiling(cap / dy)
  D <- matrix(cap^2, nx, ny)
  for (k in seq_len(nrow(z))) {
    ci <- round((z[k, 1] - gx[1]) / dx) + 1
    cj <- round((z[k, 2] - gy[1]) / dy) + 1
    ii <- max(1, ci - wx):min(nx, ci + wx)
    jj <- max(1, cj - wy):min(ny, cj + wy)
    D[ii, jj] <- pmin(D[ii, jj], outer((gx[ii] - z[k, 1])^2, (gy[jj] - z[k, 2])^2, "+"))
  }
  sqrt(D)
}
