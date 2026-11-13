## aov() models with Error() strata (class "aovlist"), analysed by strata as
## Genstat's ANOVA does, without refitting them as mixed models.
##
## The strata are those of aov(): the rows of Q'y, Q from the QR
## decomposition of the Error() model, grouped by error term (the last group
## is "Within"). In stratum i the rotated design X_i and response y_i give a
## least-squares fit whose residual mean square xi_i (df_i) estimates the
## stratum variance. The coefficients combine the strata's information,
##   b = (sum_i X_i'X_i / xi_i)^-1 sum_i X_i'y_i / xi_i,
## which for orthogonal (balanced) designs is the stratum estimate of each
## effect (and for non-orthogonal ones recovers inter-stratum information
## with ANOVA estimates of the stratum variances). Their covariance is
## V = (sum_i X_i'X_i / xi_i)^-1, so a difference of means estimated in one
## stratum gets that stratum's mean square, whatever the other strata hold
## (a negative variance component by ANOVA included, where REML would fix it
## at zero). A comparison drawing on several strata gets Satterthwaite's df,
## sum(v_i)^2 / sum(v_i^2 / df_i) with v_i the part of its variance from
## stratum i (Genstat's rule). A stratum without residual df (the grand mean,
## blocks with no error df) takes the variance of the nearest stratum below
## it that has residual df (its df too), which only matters for the means
## themselves, not for their differences.

## A carrier for the reference grid: the fixed part as an lm (its fit is not
## used for any inference), of class "pm_aovlist", holding the aovlist.
pm_aov_carrier <- function(model) {
  stopifnot(inherits(model, "aovlist"))
  Terms <- stats::terms(model)
  ind <- attr(Terms, "specials")$Error
  errorterm <- attr(Terms, "variables")[[1L + ind]]
  form <- stats::formula(Terms)
  fixed <- stats::update(form, paste(". ~ . -", paste(deparse(errorterm), collapse = "")))
  cl <- attr(model, "call")
  cl[[1L]] <- quote(stats::lm)
  cl$formula <- fixed
  cl$projections <- NULL
  fit <- tryCatch(eval(cl, environment(form)), error = function(e) NULL)
  if (is.null(fit)) {
    mf <- stats::model.frame(model)
    fit <- stats::lm(fixed, data = mf)
  }
  attr(fit, "aovlist") <- model
  attr(fit, "error_formula") <- errorterm[[2L]]
  class(fit) <- c("pm_aovlist", class(fit))
  fit
}

## The strata and the combined estimates. 'carrier' from pm_aov_carrier().
## Returns list(mp (as mymodelparm), parts (class "pm_satt", for the
## Satterthwaite df and F tests), strata (name, df, mean square)).
pm_aov_strata <- function(carrier) {
  mf <- stats::model.frame(carrier)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(carrier)
  b_lm <- stats::coef(carrier)
  est <- !is.na(b_lm)
  X <- X[, est, drop = FALSE]
  # the Error() model, as aov() builds it (Helmert contrasts)
  aovl <- attr(carrier, "aovlist")
  dat <- stats::model.frame(aovl)
  op <- options(contrasts = c("contr.helmert", "contr.poly"))
  on.exit(options(op))
  Xe <- stats::model.matrix(stats::as.formula(paste("~", paste(deparse(attr(carrier, "error_formula")), collapse = ""))),
                            dat)
  options(op)
  qe <- qr(Xe)
  rk <- qe$rank
  asg <- attr(Xe, "assign")[qe$pivot[seq_len(rk)]]
  n <- nrow(Xe)
  grp <- c(asg, rep(max(asg) + 1L, n - rk))           # stratum of each rotated row; the last is "Within"
  qtx <- qr.qty(qe, X)
  qty <- as.vector(qr.qty(qe, y))
  st <- sort(unique(grp))
  fits <- lapply(st, function(g) {
    sel <- grp == g
    Xi <- qtx[sel, , drop = FALSE]
    yi <- qty[sel]
    cols <- colSums(Xi^2) > 1e-5 * pmax(colSums(X^2), 1)
    r <- if (any(cols)) qr(Xi[, cols, drop = FALSE])$rank else 0L
    res <- if (any(cols)) stats::lm.fit(Xi[, cols, drop = FALSE], yi)$residuals else yi
    dfi <- sum(sel) - r
    list(X = Xi, y = yi, df = dfi, ms = if (dfi > 0) sum(res^2) / dfi else NA_real_)
  })
  df <- vapply(fits, `[[`, numeric(1), "df")
  ms <- vapply(fits, `[[`, numeric(1), "ms")
  # strata without residual df take the variance of the nearest stratum below
  # them with residual df (above, if none below)
  known <- which(df > 0)
  if (!length(known)) stop("no stratum of the aov model has residual degrees of freedom.", call. = FALSE)
  donor <- vapply(seq_along(fits), function(i) {
    if (df[i] > 0) return(i)
    below <- known[known > i]
    if (length(below)) min(below) else max(known[known < i])
  }, integer(1))
  xi <- ms[donor]
  info <- Reduce(`+`, lapply(seq_along(fits), function(i) crossprod(fits[[i]]$X) / xi[i]))
  score <- Reduce(`+`, lapply(seq_along(fits), function(i) crossprod(fits[[i]]$X, fits[[i]]$y) / xi[i]))
  V <- solve(info)
  b <- drop(V %*% score)
  dimnames(V) <- list(colnames(X), colnames(X))
  names(b) <- colnames(X)
  # d V / d xi_d for each stratum with residual df d (its own and the strata it lends to)
  Jac <- lapply(known, function(d) {
    ii <- which(donor == d)
    V %*% Reduce(`+`, lapply(ii, function(i) crossprod(fits[[i]]$X) / xi[i]^2)) %*% V
  })
  parts <- structure(list(vcov_beta = V, vcov_varpar = diag(2 * ms[known]^2 / df[known], nrow = length(known)),
                          Jac_list = Jac), class = "pm_satt")
  nm <- c("(Intercept)", attr(stats::terms(stats::as.formula(paste("~", paste(deparse(attr(carrier, "error_formula")),
                                                                                      collapse = "")))), "term.labels"))
  nm <- c(nm, "Within")[st + 1L]
  mp <- structure(list(coef = b, vcov = V, df = NA_real_, estimable = est), class = "mymodelparm")
  list(mp = mp, parts = parts, strata = data.frame(stratum = nm, df = df, ms = ms))
}
