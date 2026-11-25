#' Sample Variogram of the Residuals
#'
#' @description
#' The sample variogram of a model's residuals over the field: half the mean squared difference of the residuals of
#' two plots, by how far apart they are. Plots that are close together and more alike than the model allows (spatial
#' correlation, a trend across the field) show as a variogram that rises with the lag instead of staying flat. A
#' global envelope from permutations of the residuals shows the range expected without any spatial pattern.
#'
#' Two layouts:
#' \itemize{
#' \item \strong{Row x column grid} (\code{row}, \code{col}): field trials. The semivariance of plots \code{lag_row}
#'   rows and \code{lag_col} columns apart (in either diagonal direction), as ASReml and Genstat show it; the edges
#'   \code{lag_col = 0} (plots in the same column) and \code{lag_row = 0} (in the same row) are the two directions an
#'   AR1 x AR1 model describes.
#' \item \strong{Coordinates} (\code{coords}): any positions. The semivariance by class of distance (class k holds
#'   the distances in [k, k+1) x \code{step}, as Genstat's \code{FVARIOGRAM}), over all directions or by
#'   \code{directions}.
#' }
#' @param model A fitted model (any class \code{\link{residplot}} handles), or a numeric vector of data values, whose
#'   variogram is then given as it is (not scaled, no advice).
#' @param row,col Names of the variables (in \code{data}) holding the row and column positions of the plots, as whole
#'   numbers counted over the whole field (or over each trial with \code{by}), not within blocks; or vectors of them.
#' @param coords Instead of \code{row} and \code{col}: the names of two variables with the x and y coordinates (or a
#'   two-column matrix); with \code{lonlat = TRUE} longitude and latitude in degrees.
#' @param by Name of a variable (or a vector) that splits the plots into trials: one variogram per trial, never pooled.
#' @param data The data: by default those the model was fitted to (the rows used); needed for a numeric
#'   \code{model}.
#' @param max_lag Grid: the largest row and column lags (one number or two); default half the number of rows and of
#'   columns (longer lags have few pairs and add noise to the test). Coordinates: the largest distance; default half the largest distance between two points.
#' @param step Coordinates: the width of the distance classes; default \code{max_lag / 15}.
#' @param directions Coordinates: directions in degrees anticlockwise from the x axis (east), e.g.
#'   \code{c(0, 45, 90, 135)}; default all directions together.
#' @param segment Coordinates: the angle (degrees) each direction covers, centred on it; default 45.
#' @param lonlat \code{TRUE}: \code{coords} are longitude and latitude (degrees); distances in km (haversine).
#' @param robust \code{TRUE}: Cressie and Hawkins' estimator, less affected by a few large differences.
#' @param min_pairs Lags with fewer pairs are left out (default 20, as Genstat).
#' @param nsim Number of data sets simulated from the fitted model for the envelope (permutations within each trial
#'   for a numeric \code{model}); 0 for none.
#' @param conf_level Coverage of the global envelope.
#' @param seed Seed of the permutations (the user's random number stream is left as it was).
#' @param ncore Number of cores for the permutations.
#' @param plot \code{TRUE}: draw the variogram (\code{autoplot()}).
#' @details
#' Residuals are those of \code{\link{residplot}} (one place): conditional residuals scaled by the residual
#' standard deviation for \code{lmer} (the random effects, e.g. blocks, removed), normalized for \code{lme} and
#' \code{gls} (their correlation structure removed: a correct spatial model leaves a flat variogram), randomized
#' quantile residuals for Poisson, binomial and negative binomial models. They are scaled by their standard deviation
#' within each trial, so the semivariance is relative to the residual variance: about 1 at every lag when there is no
#' spatial pattern (Genstat's "0-1 scale").
#'
#' The semivariance is \eqn{\gamma(h) = \sum (e_i - e_j)^2 / (2 N(h))} over the \eqn{N(h)} pairs at lag \eqn{h}
#' (Matheron); \code{robust = TRUE} gives \eqn{(\sum |e_i - e_j|^{1/2} / N)^4 / (2 (0.457 + 0.494 / N))} (Cressie and
#' Hawkins 1980). The envelope: \code{nsim} data sets are simulated from the fitted model (independent errors, no
#' spatial pattern beyond the model), the model refitted and the variogram of its residuals computed; for \code{lm}
#' and \code{aov} exactly, without refits. Residuals are correlated through the fitted effects (plots in one block
#' sum to zero), so permuting them would give a wrong reference (in a simulated trial with blocks, 15\% false alarms
#' at the 5\% level); the simulated residuals carry the same correlation. The envelope is the global scaled maximum
#' absolute difference envelope of Myllymäki et al. (2017): the observed variogram lies wholly inside it with
#' probability \code{conf_level} when the model is right, however many lags are shown; \code{p.value} is the
#' matching global test. Classes without \code{simulate()} and a refit (\code{gam}, \code{semireg}, \code{aov} with
#' \code{Error()}) get the variogram without an envelope. The test informs; the remedy is a model choice.
#'
#' A row x column position may hold one plot per trial: positions that repeat stop with a message (rows or columns
#' numbered within blocks give wrong lags).
#' @return A data frame of class \code{"pm_variogram"}: \code{trial} (with \code{by}), \code{lag_row} and
#'   \code{lag_col} (grid) or \code{direction} and \code{lag} (the mean distance in the class), \code{gamma},
#'   \code{n_pairs}, \code{lower} and \code{upper} (the envelope), \code{outside}; attributes \code{test} (the global
#'   test per trial) and advice (\code{\link{advice}}). \code{autoplot()} draws it.
#' @references
#' Cressie, N. and Hawkins, D. M. (1980). Robust estimation of the variogram. \emph{Mathematical Geology}, 12,
#' 115-125.
#'
#' Gilmour, A. R., Cullis, B. R. and Verbyla, A. P. (1997). Accounting for natural and extraneous variation in the
#' analysis of field experiments. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, 2,
#' 269-293.
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial
#' processes. \emph{Journal of the Royal Statistical Society B}, 79, 381-404.
#'
#' Webster, R. and Oliver, M. A. (1990). \emph{Statistical Methods in Soil and Land Resource Survey}. Oxford
#' University Press.
#' @seealso \code{\link{residplot}}, \code{\link{pm_check}} (with \code{row}, \code{col} or \code{coords})
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' # a field trial in 6 rows x 6 columns, 3 blocks of 2 rows, with a smooth trend over the field
#' set.seed(1)
#' d <- expand.grid(row = 1:6, col = 1:6)
#' d$block <- factor((d$row + 1) %/% 2)
#' d$variety <- factor(sample(rep(1:6, 6)))
#' d$yield <- 10 + as.numeric(d$variety) + sin(d$row / 2) + cos(d$col / 2) + rnorm(36, sd = 0.3)
#' fm <- lm(yield ~ block + variety, data = d)
#' v <- residvariogram(fm, row = "row", col = "col", nsim = 19, plot = FALSE)
#' v
#' ggplot2::autoplot(v)
#' @export
residvariogram <- function(model, row = NULL, col = NULL, coords = NULL, by = NULL, data = NULL, max_lag = NULL,
                           step = NULL, directions = NULL, segment = 45, lonlat = FALSE, robust = FALSE,
                           min_pairs = 20, nsim = 99, conf_level = 0.95, seed = 2026, ncore = 1, plot = TRUE) {
  grid <- !is.null(row) || !is.null(col)
  if (grid == !is.null(coords)) stop("Give either 'row' and 'col' (a row x column grid) or 'coords'.", call. = FALSE)
  if (grid && (is.null(row) || is.null(col))) stop("Give both 'row' and 'col'.", call. = FALSE)
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim < 0) stop("'nsim' must be 0 or a positive number.", call. = FALSE)
  raw <- is.numeric(model) && is.null(dim(model))
  if (raw) {
    e <- as.numeric(model)
    d <- data
    rlab <- "Values"
  } else {
    d <- if (is.null(data)) pm_fit_rows(model) else as.data.frame(data)
    model <- pm_fit(model)
    if (inherits(model, "lmerModLmerTest")) model <- methods::as(model, "lmerMod")
    rd <- suppressWarnings(suppressMessages(pm_resid_data(model)))
    e <- rd$resid
    rlab <- attr(rd, "rlab")
    if (is.data.frame(d) && nrow(d) != length(e)) {
      stop("The data have ", nrow(d), " rows but the model has ", length(e), " residuals; give the rows used in the ",
           "fit as 'data'.", call. = FALSE)
    }
  }
  n <- length(e)
  get <- function(v, what) pm_vario_var(v, d, n, what)
  trial <- if (is.null(by)) factor(rep("all", n)) else factor(get(by, "by"))
  ok <- is.finite(e) & !is.na(trial)
  if (grid) {
    pos <- cbind(pm_vario_position(get(row, "row"), "row"), pm_vario_position(get(col, "col"), "col"))
  } else {
    pos <- if (is.matrix(coords) || is.data.frame(coords)) as.matrix(coords) else
      do.call(cbind, lapply(coords, function(v) get(v, "coords")))
    if (ncol(pos) != 2 || nrow(pos) != n || !is.numeric(pos)) {
      stop("'coords' must name two numeric variables (x and y), or be a two-column matrix with one row per ",
           "observation.", call. = FALSE)
    }
  }
  ok <- ok & stats::complete.cases(pos)
  if (sum(ok) < 3) stop("Fewer than 3 plots with a residual and a position.", call. = FALSE)
  levs <- levels(droplevels(trial[ok]))

  parts <- lapply(levs, function(t) {
    k <- which(ok & trial == t)
    ek <- e[k]
    if (!raw) ek <- ek / stats::sd(ek)
    lab <- if (is.null(by)) NULL else t
    setup <- if (grid) pm_vario_grid_setup(pos[k, , drop = FALSE], max_lag, min_pairs, lab) else
      pm_vario_dist_setup(pos[k, , drop = FALSE], max_lag, step, directions, segment, lonlat, min_pairs, lab)
    list(e = ek, setup = setup, trial = t)
  })
  stat <- function(p, ek) if (grid) pm_vario_grid_gamma(p$setup, ek, robust) else pm_vario_dist_gamma(p$setup, ek, robust)
  # the reference: residuals of data simulated from the fitted model (values: permutations)
  sims <- if (nsim > 0) pm_vario_sims(if (raw) NULL else model, e, trial, ok, nsim, seed, ncore) else NULL
  if (nsim > 0 && is.null(sims)) {
    message("No envelope: residvariogram() simulates from lm, aov, glm, lmer, glmer, glmmTMB, lme and gls fits.")
    nsim <- 0
  }

  out <- lapply(seq_along(parts), function(i) {
    p <- parts[[i]]
    g <- stat(p, p$e)
    tab <- p$setup$table
    tab$gamma <- g
    env <- if (nsim > 0) pm_vario_envelope(p, g, stat, sims[, which(ok & trial == p$trial), drop = FALSE], raw,
                                           conf_level) else NULL
    tab$lower <- if (is.null(env)) NA_real_ else env$lower
    tab$upper <- if (is.null(env)) NA_real_ else env$upper
    tab$outside <- if (is.null(env)) NA else tab$gamma < tab$lower | tab$gamma > tab$upper
    if (!is.null(by)) tab <- cbind(trial = p$trial, tab)
    list(tab = tab, p = if (is.null(env)) NA_real_ else env$p, missing = p$setup$missing)
  })
  res <- do.call(rbind, lapply(out, `[[`, "tab"))
  rownames(res) <- NULL
  test <- data.frame(trial = levs, n = vapply(parts, function(p) length(p$e), integer(1)),
                     p.value = vapply(out, `[[`, numeric(1), "p"),
                     missing = vapply(out, `[[`, numeric(1), "missing"), stringsAsFactors = FALSE)
  res <- structure(res, class = c("pm_variogram", "data.frame"), layout = if (grid) "grid" else "distance",
                   test = test, nsim = if (is.null(sims)) 0L else nrow(sims), conf_level = conf_level, robust = robust, rlab = rlab, raw = raw,
                   lonlat = lonlat, by = !is.null(by))
  if (!raw) res <- pm_add_advice(res, pm_variogram_advice(res))
  if (plot) print(autoplot(res))
  res
}

## A variable of the data by name, or a vector given directly (one value per observation).
pm_vario_var <- function(v, d, n, what) {
  if (is.character(v) && length(v) == 1) {
    if (!is.data.frame(d) || !v %in% names(d)) {
      stop("'", what, "': there is no variable \"", v, "\" in the data; give 'data'.", call. = FALSE)
    }
    return(d[[v]])
  }
  if (length(v) != n) stop("'", what, "' must be a variable name or have one value per observation (", n, ").",
                           call. = FALSE)
  v
}

## Row or column positions as whole numbers (factors with numeric labels are read as their numbers).
pm_vario_position <- function(x, what) {
  if (is.factor(x) || is.character(x)) {
    v <- suppressWarnings(as.numeric(as.character(x)))
    if (anyNA(v[!is.na(x)])) {
      stop("'", what, "' must hold numbers (row or column positions); \"", x[!is.na(x) & is.na(v)][1],
           "\" is not one.", call. = FALSE)
    }
    x <- v
  }
  if (!is.numeric(x) || any(abs(x - round(x)) > 1e-8, na.rm = TRUE)) {
    stop("'", what, "' must hold whole numbers (row or column positions); for other positions use 'coords'.",
         call. = FALSE)
  }
  round(x)
}

## Grid lags of one trial: the residuals in an nr x nc matrix (NA where no
## plot); for each absolute lag (dr, dc) the pairs are the shifted matrices,
## in both diagonal directions when dr > 0 and dc > 0. The pairs counts depend
## only on which positions hold a plot, so the lags kept (at least min_pairs
## pairs) are found once and the permutations reuse them.
pm_vario_grid_setup <- function(pos, max_lag, min_pairs, trial = NULL) {
  r <- pos[, 1] - min(pos[, 1]) + 1
  c <- pos[, 2] - min(pos[, 2]) + 1
  key <- paste(r, c)
  if (anyDuplicated(key)) {
    dup <- pos[duplicated(key), , drop = FALSE][1, ]
    stop("Two or more plots share row ", dup[1], " and column ", dup[2], if (!is.null(trial)) paste0(" in trial ", trial),
         ": 'row' and 'col' must number the rows and columns of the whole field",
         if (is.null(trial)) " (or of each trial, with 'by')", ", not within blocks. If they count within blocks, ",
         "make field positions first (e.g. row + (block - 1) * rows per block).", call. = FALSE)
  }
  nr <- max(r); nc <- max(c)
  ml <- if (is.null(max_lag)) pmax(1, c(nr, nc) %/% 2) else rep_len(as.integer(max_lag), 2)
  ml <- pmin(ml, c(nr - 1, nc - 1))
  idx <- cbind(r, c)
  occ <- matrix(0, nr, nc)
  occ[idx] <- 1
  lags <- expand.grid(lag_row = 0:ml[1], lag_col = 0:ml[2])[-1, ]
  cnt <- pm_vario_grid_sums(occ, lags, function(a, b) a * b)
  keep <- cnt >= min_pairs
  lags <- lags[keep, , drop = FALSE]
  rownames(lags) <- NULL
  list(table = cbind(lags, n_pairs = cnt[keep]), idx = idx, nr = nr, nc = nc,
       missing = 1 - nrow(pos) / (nr * nc))
}

## For each lag, the sum of f(a, b) over the pairs of grid cells that lag
## apart (a and b the two cells; NA cells give NA and are skipped).
pm_vario_grid_sums <- function(M, lags, f) {
  nr <- nrow(M); nc <- ncol(M)
  vapply(seq_len(nrow(lags)), function(i) {
    dr <- lags$lag_row[i]; dc <- lags$lag_col[i]
    if (dr >= nr || dc >= nc) return(0)
    a <- M[seq_len(nr - dr), seq_len(nc - dc), drop = FALSE]
    s <- sum(f(a, M[dr + seq_len(nr - dr), dc + seq_len(nc - dc), drop = FALSE]), na.rm = TRUE)
    if (dr > 0 && dc > 0) {
      s <- s + sum(f(M[seq_len(nr - dr), dc + seq_len(nc - dc), drop = FALSE],
                     M[dr + seq_len(nr - dr), seq_len(nc - dc), drop = FALSE]), na.rm = TRUE)
    }
    s
  }, numeric(1))
}

## The semivariances at the lags of a grid setup, for residuals e.
pm_vario_grid_gamma <- function(setup, e, robust = FALSE) {
  E <- matrix(NA_real_, setup$nr, setup$nc)
  E[setup$idx] <- e
  lags <- setup$table
  N <- lags$n_pairs
  if (!robust) return(pm_vario_grid_sums(E, lags, function(a, b) (a - b)^2) / (2 * N))
  m <- pm_vario_grid_sums(E, lags, function(a, b) sqrt(abs(a - b))) / N
  m^4 / (2 * (0.457 + 0.494 / N))
}

## Distance classes of one set of points: the pairs within max_lag (and the
## directions), each with its class; kept once for the permutations.
pm_vario_dist_setup <- function(pos, max_lag, step, directions, segment, lonlat, min_pairs, trial = NULL) {
  n <- nrow(pos)
  if (n * (n - 1) / 2 > 2e7) stop("Too many points (", n, ") for all pairs; use a subset or 'by'.", call. = FALSE)
  pr <- which(upper.tri(matrix(TRUE, n, n)), arr.ind = TRUE)
  i <- pr[, 1]; j <- pr[, 2]
  if (lonlat) {
    rad <- pi / 180
    lon <- pos[, 1] * rad; lat <- pos[, 2] * rad
    a <- sin((lat[j] - lat[i]) / 2)^2 + cos(lat[i]) * cos(lat[j]) * sin((lon[j] - lon[i]) / 2)^2
    h <- 2 * 6371.0088 * asin(pmin(1, sqrt(a)))
    # direction on the local plane: east and north displacements
    dx <- (lon[j] - lon[i]) * cos((lat[i] + lat[j]) / 2); dy <- lat[j] - lat[i]
  } else {
    dx <- pos[j, 1] - pos[i, 1]; dy <- pos[j, 2] - pos[i, 2]
    h <- sqrt(dx^2 + dy^2)
  }
  if (is.null(max_lag)) max_lag <- max(h) / 2
  if (is.null(step)) step <- max_lag / 15
  if (!is.numeric(step) || step <= 0 || !is.numeric(max_lag) || max_lag <= 0) {
    stop("'step' and 'max_lag' must be positive.", call. = FALSE)
  }
  k <- floor(h / step)                                     # class k: [k, k+1) x step (Genstat FVARIOGRAM)
  use <- h > 0 & h < max_lag
  ang <- (atan2(dy, dx) * 180 / pi) %% 180
  dirs <- if (is.null(directions)) NA_real_ else directions %% 180
  sets <- lapply(dirs, function(dir) {
    sel <- use
    if (!is.na(dir)) sel <- sel & abs(((ang - dir + 90) %% 180) - 90) <= segment / 2
    cls <- k[sel]
    g <- factor(cls, levels = sort(unique(cls)))
    tab <- data.frame(direction = dir, class = as.integer(levels(g)),
                      lag = as.vector(tapply(h[sel], g, mean)), n_pairs = as.vector(table(g)))
    keep <- tab$n_pairs >= min_pairs
    list(tab = tab[keep, , drop = FALSE], i = i[sel], j = j[sel], g = g, keep = keep)
  })
  tab <- do.call(rbind, lapply(sets, `[[`, "tab"))
  if (is.null(directions)) tab$direction <- NULL
  tab$class <- NULL
  rownames(tab) <- NULL
  list(table = tab, sets = sets, missing = 0)
}

## The semivariances of the distance classes of a setup, for values e.
pm_vario_dist_gamma <- function(setup, e, robust = FALSE) {
  unlist(lapply(setup$sets, function(s) {
    dz <- e[s$i] - e[s$j]
    if (!robust) return((tapply(dz^2, s$g, mean) / 2)[s$keep])
    N <- as.vector(table(s$g))
    m <- as.vector(tapply(sqrt(abs(dz)), s$g, mean))
    (m^4 / (2 * (0.457 + 0.494 / N)))[s$keep]
  }), use.names = FALSE)
}

## The reference of the envelope: an nsim x n matrix of residuals of data
## simulated from the fitted model (no spatial pattern beyond the model),
## so that the correlation the fit itself gives the residuals (e.g. within
## blocks) is in the reference too. lm/aov: exactly, D (I - H) z with D the
## scaling of rstandard() (the scale is removed per trial later). Other
## classes: simulate() and refit (R/boot-parts.R), residuals as residplot().
## Values (model NULL): permutations within each trial. NULL: no reference.
pm_vario_sims <- function(model, e, trial, ok, nsim, seed, ncore) {
  n <- length(e)
  nsim <- as.integer(nsim)
  if (is.null(model)) {
    return(pm_with_seed(seed, t(vapply(seq_len(nsim), function(b) {
      out <- e
      for (t in levels(trial)) {
        k <- which(ok & trial == t)
        out[k] <- e[k][sample.int(length(k))]
      }
      out
    }, numeric(n)))))
  }
  if (inherits(model, "lm") && !inherits(model, c("glm", "gam", "mlm")) && !is.null(model$qr)) {
    h <- stats::hatvalues(model)
    Z <- pm_with_seed(seed, matrix(stats::rnorm(length(h) * nsim), length(h)))
    R <- base::qr.resid(model$qr, Z) / sqrt(pmax(1 - h, 0))
    R[h > 1 - 1e-10, ] <- NA
    return(if (nrow(R) == n) t(R) else NULL)
  }
  if (!inherits(model, c("merMod", "lme", "gls", "glmmTMB", "glm")) || inherits(model, "gam")) return(NULL)
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  f <- stats::formula(model)
  new_lhs <- function(f) {
    f[[2]] <- quote(.pm_y)
    f
  }
  bp <- tryCatch(pm_boot_fit(model, f, new_lhs, quiet), error = function(err) NULL)  # (its "null" fit: the model)
  if (is.null(bp)) return(NULL)
  ysim <- pm_with_seed(seed, bp$simulate(nsim))
  res <- pm_perm_map(ysim, function(y) {
    r <- quiet(pm_resid_data(bp$refit(y)))$resid
    if (length(r) == n) r else NULL
  }, ncore = ncore, packages = unique(c("predictmeans", bp$packages)))
  res <- pm_perm_drop(res)$res
  if (!length(res)) return(NULL)
  do.call(rbind, res)
}

## Global envelope (scaled maximum absolute difference, Myllymäki et al.
## 2017) of one trial from its simulated residuals S (nsim x n): the
## variogram stays inside it at every lag with probability conf_level when
## the model is right; p is the matching global test.
pm_vario_envelope <- function(p, g, stat, S, raw, conf_level) {
  G <- vapply(seq_len(nrow(S)), function(b) {
    eb <- S[b, ]
    if (!raw) eb <- eb / stats::sd(eb, na.rm = TRUE)
    stat(p, eb)
  }, numeric(length(g)))
  G <- matrix(G, nrow = length(g))
  if (!length(g)) return(NULL)
  # centre and spread from all B + 1 curves, so the observed one is exchangeable with the simulated (exact p)
  A <- cbind(g, G)
  mu <- rowMeans(A, na.rm = TRUE)
  s <- apply(A, 1, stats::sd, na.rm = TRUE)
  s[!is.finite(s) | s <= 0] <- NA
  Tsim <- apply(abs(G - mu) / s, 2, max, na.rm = TRUE)
  Tobs <- max(abs(g - mu) / s, na.rm = TRUE)
  B <- length(Tsim)
  # the k-th largest of the simulated maxima, k = floor(alpha (B + 1))
  k <- max(1L, floor((1 - conf_level) * (B + 1)))
  u <- sort(Tsim, decreasing = TRUE)[k]
  s[is.na(s)] <- 0
  list(lower = pmax(mu - u * s, 0), upper = mu + u * s, p = (1 + sum(Tsim >= Tobs)) / (B + 1))
}

## print() of residvariogram(): the global test, the short lags, the advice.
#' @exportS3Method print pm_variogram
print.pm_variogram <- function(x, digits = 3, ...) {
  test <- attr(x, "test")
  if (is.null(test)) {
    print(structure(x, class = "data.frame"), digits = digits, row.names = FALSE)
    return(invisible(x))
  }
  grid <- attr(x, "layout") == "grid"
  cat("Sample variogram of ", if (isTRUE(attr(x, "raw"))) "the values" else
    tolower(attr(x, "rlab")), if (grid) " on the row x column grid" else " by distance",
    if (isTRUE(attr(x, "robust"))) " (Cressie-Hawkins)" else "", "\n", sep = "")
  if (!isTRUE(attr(x, "raw"))) cat("Semivariance relative to the residual variance: about 1 without spatial pattern.\n")
  if (attr(x, "nsim") > 0) {
    pv <- format.pval(test$p.value, digits = 2, eps = 1e-4)
    cat("Global envelope test (", format(100 * attr(x, "conf_level")), "%, ", attr(x, "nsim"),
        if (isTRUE(attr(x, "raw"))) " permutations" else " simulations from the model", "): ",
        if (isTRUE(attr(x, "by"))) paste0(test$trial, " p = ", pv, collapse = "; ") else paste0("p = ", pv),
        "\n", sep = "")
  }
  y <- structure(x, class = "data.frame")
  for (a in c("layout", "test", "nsim", "conf_level", "robust", "rlab", "raw", "lonlat", "by", "advice")) attr(y, a) <- NULL
  if (grid) {
    edge <- y$lag_row == 0 | y$lag_col == 0
    cat("\nLags along rows and columns (all", nrow(y), "lags: as.data.frame(x); the surface: autoplot(x, type =",
        "\"surface\")):\n")
    y <- y[edge, , drop = FALSE]
    y <- y[order(if (is.null(y$trial)) rep(0, nrow(y)) else y$trial, y$lag_col > 0, y$lag_row + y$lag_col), , drop = FALSE]
  } else cat("\n")
  print(y, digits = digits, row.names = FALSE)
  pm_print_advice(x)
  invisible(x)
}

#' Plot a Residual Variogram
#'
#' @description
#' \code{autoplot()} of \code{\link{residvariogram}}: the semivariance against the lag with the global envelope
#' shaded and the lags outside it marked in red. Grid layouts: \code{type = "edges"} (default) the lags along a
#' column and along a row; \code{"surface"} every row x column lag as tiles.
#' @param object A result of \code{residvariogram()}.
#' @param type \code{"edges"} or \code{"surface"} (grid layouts only).
#' @param ... Not used.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{residvariogram}}
#' @exportS3Method ggplot2::autoplot pm_variogram
autoplot.pm_variogram <- function(object, type = c("edges", "surface"), ...) {
  type <- match.arg(type)
  st <- pm_curve_style_defaults
  d <- structure(object, class = "data.frame")
  grid <- attr(object, "layout") == "grid"
  by <- isTRUE(attr(object, "by"))
  ylab <- if (isTRUE(attr(object, "raw"))) "Semivariance" else "Semivariance (relative to residual variance)"
  ref <- if (isTRUE(attr(object, "raw"))) NULL else ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                                                                         colour = "grey50")
  th <- ggplot2::theme_bw(base_size = 10)
  if (grid && type == "surface") {
    p <- ggplot2::ggplot(d, ggplot2::aes(.data$lag_col, .data$lag_row, fill = .data$gamma)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::scale_fill_gradient(low = "#F7FBFF", high = st$line, name = "Semi-\nvariance") +
      ggplot2::scale_y_reverse() +
      ggplot2::labs(x = "Column lag", y = "Row lag") + th
    if (any(d$outside, na.rm = TRUE)) {
      p <- p + ggplot2::geom_point(data = d[which(d$outside), , drop = FALSE], colour = st$ref, size = 1.5,
                                   show.legend = FALSE)
    }
    if (by) p <- p + ggplot2::facet_wrap(~trial)
    return(p)
  }
  if (grid) {
    d <- d[d$lag_row == 0 | d$lag_col == 0, , drop = FALSE]
    d$along <- factor(ifelse(d$lag_col == 0, "Same column (row lag)", "Same row (column lag)"),
                      levels = c("Same column (row lag)", "Same row (column lag)"))
    d$lag <- d$lag_row + d$lag_col
    panel <- if (by) ggplot2::facet_grid(trial ~ along) else ggplot2::facet_wrap(~along)
    xlab <- "Lag (rows or columns)"
  } else {
    panel <- if (by && !is.null(d$direction)) ggplot2::facet_grid(trial ~ direction, labeller = ggplot2::label_both)
    else if (by) ggplot2::facet_wrap(~trial) else if (!is.null(d$direction))
      ggplot2::facet_wrap(~direction, labeller = ggplot2::label_both)
    xlab <- if (isTRUE(attr(object, "lonlat"))) "Distance (km)" else "Distance"
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$lag, .data$gamma)) + ref
  if (!all(is.na(d$lower))) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), fill = st$band,
                                  alpha = st$alpha)
  }
  p <- p + ggplot2::geom_line(colour = st$line) + ggplot2::geom_point(colour = st$line, size = 1.3)
  if (any(d$outside, na.rm = TRUE)) {
    p <- p + ggplot2::geom_point(data = d[which(d$outside), , drop = FALSE], colour = st$ref, size = 2)
  }
  p + panel + ggplot2::labs(x = xlab, y = ylab) + th
}
