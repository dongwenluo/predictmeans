## Overall means (1 - pi) mu of glmmTMB fits with zero inflation, for
## predictmeans(type = "response").
##
## For each cell of the reference grid (rg, from Kmatrix(refgrid = TRUE)) the
## mean of the count part is mu = exp(eta) (log link) and the probability of
## an extra zero pi = plogis(zeta) (eta, zeta: the fixed parts of the
## conditional and zero-inflation models, random effects at zero); the
## overall mean is (1 - pi) mu, averaged over the cells with the weights of
## the means (as average = "response"). Its gradient with respect to the
## fixed coefficients of both parts: (1 - pi) mu x for the conditional ones,
## -mu pi (1 - pi) x_zi for the zero-inflation ones; their joint covariance
## from vcov(full = TRUE).

## Rows 'rows' of rg$W: returns list(m, J (gradient), V (covariance of the
## conditional and zero-inflation fixed coefficients)).
pm_zi_overall <- function(model, rg, rows) {
  fx <- glmmTMB::fixef(model)
  bc <- fx$cond
  bz <- fx$zi
  nc <- length(bc)
  nz <- length(bz)
  Vf <- as.matrix(stats::vcov(model, full = TRUE))
  V <- Vf[seq_len(nc + nz), seq_len(nc + nz), drop = FALSE]
  grid <- rg$grid
  X <- as.matrix(rg$basis(grid))[, names(bc), drop = FALSE]
  fz <- stats::formula(model, fixed.only = TRUE, component = "zi")
  Xz <- tryCatch({
    tz <- stats::delete.response(stats::terms(fz))
    stats::model.matrix(tz, stats::model.frame(tz, grid, na.action = stats::na.pass))
  }, error = function(e) NULL)
  if (is.null(Xz) || !all(names(bz) %in% colnames(Xz))) {
    stop("type = \"response\": the zero-inflation formula (", paste(deparse(fz), collapse = ""), ") uses variables or ",
         "terms that are not in the grid of the means; give them in 'modelterm' or 'covariate'.", call. = FALSE)
  }
  Xz <- Xz[, names(bz), drop = FALSE]
  mu <- exp(as.vector(X %*% bc))
  p0 <- stats::plogis(as.vector(Xz %*% bz))
  W <- as.matrix(rg$W[rows, , drop = FALSE])
  m <- as.vector(W %*% ((1 - p0) * mu))
  J <- cbind(W %*% (((1 - p0) * mu) * X), W %*% ((-mu * p0 * (1 - p0)) * Xz))
  list(m = m, J = J, V = V)
}

## Replaces the estimates of pm_estimate() by those of the log overall means:
## the K matrix becomes the identity over the means, the coefficients the log
## means and their covariance J V J' (on the log scale), so the comparisons,
## SEDs, letters and the back-transformation (exp, the inverse link) of the
## later stages work unchanged; the df of a contrast L of the means are those
## of L K for the conditional model (its df rule, by term).
pm_zi_response_est <- function(opt, est) {
  K0 <- est$K
  KK <- est$KK
  rows <- match(rownames(K0), rownames(KK))
  o <- pm_zi_overall(opt$model, est$ctr.matrix$refgrid, rows)
  nm <- rownames(K0)
  Jl <- o$J / o$m
  C <- Jl %*% o$V %*% t(Jl)
  dimnames(C) <- list(nm, nm)
  Kn <- diag(length(nm))
  dimnames(Kn) <- list(nm, nm)
  mp <- est$mp
  mp$coef <- stats::setNames(log(o$m), nm)
  mp$vcov <- C
  mp$estimable <- rep(TRUE, length(nm))
  inf <- est$inf
  df0 <- inf$df
  inf$df <- function(L, terms) df0(as.matrix(L) %*% K0, terms)
  inf$mp <- mp
  fct <- est$ctr.matrix$fctnames[rows, , drop = FALSE]
  e <- pm_estimate_means(Kn, fct, est$mdf, opt$vars, mp, opt$ndecimal)
  est$ctr.matrix$K <- Kn
  est$ctr.matrix$fctnames <- fct
  utils::modifyList(est, list(ctr.matrix = est$ctr.matrix, KK = Kn, inf = inf, mp = mp, vcovm = C, K = e$K,
                              label = e$label, pm = e$pm, ses = e$ses, ndf = e$ndf, mt = e$mt, bkmt = e$mt,
                              mean.table = e$mean.table, se.table = e$se.table))
}
