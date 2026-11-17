## Overall means (1 - pi) mu of glmmTMB fits with zero inflation, for
## predictmeans(mean_type = "overall").
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
    stop("mean_type = \"overall\": the zero-inflation formula (", paste(deparse(fz), collapse = ""), ") uses variables or ",
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
pm_zi_overall_est <- function(opt, est) {
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

## The covariance of the fixed and random coefficients of the conditional
## model of a glmmTMB fit, (beta, b) (the first nb random effects), given the
## other parameters (variances, zero-inflation and dispersion coefficients):
## the inverse of their block of TMB's joint precision. With zi = TRUE the
## fixed zero-inflation coefficients are included (after b), for the overall
## means of semipred(mean_type = "overall").
## full = TRUE: the block of the inverse of the whole joint precision instead
## (the uncertainty of the variance parameters included: unconditional).
pm_tmb_cond_cov <- function(m, nb, zi = FALSE, full = FALSE) {
  pm_need("TMB", "the covariance of zero-inflated glmmTMB fits")
  P <- TMB::sdreport(m$obj, getJointPrecision = TRUE)$jointPrecision
  nm <- rownames(P)
  ib <- which(nm == "beta")
  iu <- which(nm == "b")[seq_len(nb)]
  iz <- if (zi) which(nm == "betazi") else integer(0)
  idx <- c(ib, iu, iz)
  if (full) as.matrix(solve(P))[idx, idx, drop = FALSE] else as.matrix(solve(P[idx, idx]))
}

## The fixed design of the zero-inflation model at the rows of 'grid' (a data
## frame with the variables of the grid), columns 'names_z'.
pm_zi_X <- function(model, grid, names_z) {
  fz <- stats::formula(model, fixed.only = TRUE, component = "zi")
  Xz <- tryCatch({
    tz <- stats::delete.response(stats::terms(fz))
    stats::model.matrix(tz, stats::model.frame(tz, as.data.frame(grid), na.action = stats::na.pass))
  }, error = function(e) NULL)
  if (is.null(Xz) || !all(names_z %in% colnames(Xz))) {
    stop("mean_type = \"overall\": the zero-inflation formula (", paste(deparse(fz), collapse = ""), ") uses ",
         "variables or terms that are not in the grid of the curves; give them in 'modelterm' or 'covariate'.",
         call. = FALSE)
  }
  Xz[, names_z, drop = FALSE]
}

## Curves of overall means (semipred(), covariatemeans()): rows K of the
## count model (coefficients b), rows Kz of the zero-inflation model
## (coefficients bz), V the joint covariance of (b, bz). log m = eta +
## log(1 - pi), pi = plogis(zeta); returned as a linear predictor with an
## extra fixed column: K' = [K, -pi Kz, c] with c = log(1 - pi) + pi zeta and
## coefficients (b, bz, 1), so K' coef = log m exactly and K' V' K' (V' = V
## and 0 for the constant) its delta-method variance; the inverse link exp
## gives m and its derivative the SE of m, differences the log ratios of
## overall means. ncond: the columns of the count model (for the df).
pm_zi_augment <- function(K, b, Kz, bz, V) {
  K <- as.matrix(K)
  zeta <- as.vector(Kz %*% bz)
  p0 <- stats::plogis(zeta)
  Ka <- cbind(K, -p0 * Kz, log1p(-p0) + p0 * zeta)
  V <- as.matrix(V)
  Va <- rbind(cbind(V, 0), 0)
  nm <- c(names(b), paste0("zi~", names(bz)), "(log(1 - pi))")
  dimnames(Va) <- list(nm, nm)
  colnames(Ka) <- nm
  list(K = Ka, coef = stats::setNames(c(b, bz, 1), nm), vcov = Va, ncond = ncol(K))
}

## The checks of mean_type = "overall" (one place): returns the type to use.
pm_zi_check <- function(model, mean_type, deriv = 0) {
  if (!identical(mean_type, "overall")) return("conditional")
  if (is.null(pm_zi_component(model))) {
    message("mean_type = \"overall\" applies to glmmTMB models with zero inflation; this model has none, so the ",
            "means are the usual ones.")
    return("conditional")
  }
  m <- if (inherits(model, "semireg")) model$semer else model
  fam <- stats::family(m)
  if (!identical(fam$link, "log") || grepl("truncated", fam$family)) {
    stop("mean_type = \"overall\" is available for zero-inflated count models with a log link (poisson, nbinom1, ",
         "nbinom2, ...), not ", fam$family, " (", fam$link, ").", call. = FALSE)
  }
  if (deriv > 0) stop("mean_type = \"overall\" is for the curves, not their derivatives.", call. = FALSE)
  "overall"
}

## The zero-inflation message of the curve functions, by mean_type.
pm_zi_message <- function(model, mean_type) {
  zi <- pm_zi_component(model)
  if (is.null(zi)) return(invisible(NULL))
  message(if (identical(mean_type, "overall"))
    "The model has zero inflation: the curves are of the overall means (1 - pi) mu (mean_type = \"overall\")." else
    paste(attr(zi, "note"), "mean_type = \"overall\" gives the overall means."))
}
