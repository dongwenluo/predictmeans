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
  # missing values of an orthogonal design: estimated as Genstat's ANOVA does
  comp <- tryCatch(pm_aov_complete(carrier), error = function(e) NULL)
  if (!is.null(comp)) {
    res <- tryCatch(pm_aov_missing(carrier, comp, X, y), error = function(e) NULL)
    if (!is.null(res)) return(res)
  }
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
  if (!pm_aov_orthogonal(qe, grp, carrier, dat)) {
    message("aov() with Error() on a design that is not orthogonal: the strata are combined with the ANOVA ",
            "estimates of their variances. REML (e.g. lmer() with the Error() terms as random effects) uses the ",
            "information in all strata and is the better analysis.")
  }
  list(mp = mp, parts = parts, strata = data.frame(stratum = nm, df = df, ms = ms))
}

## ---- Missing values (Genstat's ANOVA) ----
##
## An orthogonal design with some plots missing: Genstat's ANOVA estimates the
## missing values so that the residual sum of squares of the lowest stratum is
## smallest (what that leaves undetermined, e.g. a whole plot missing, by the
## next stratum up, and so on), takes the residual df of each stratum down by
## the values it determines, and analyses the completed design by strata.
## The estimates are linear in the observed values (y_missing = M y_observed),
## so the coefficients are b = B y_observed and their covariance under the
## strata model is exact: B (sum_i xi_i P_i[obs, obs]) B', P_i the projector of
## stratum i; Genstat's SEDs instead treat the design as complete ("not
## adjusted for missing values").

## The units of the complete design: list(full (data frame of the design
## variables, every unit), obs (row of the model frame for each unit, NA when
## missing)), or NULL when nothing is missing or the design cannot be completed.
pm_aov_complete <- function(carrier) {
  aovl <- attr(carrier, "aovlist")
  mf <- stats::model.frame(aovl)
  ef <- attr(carrier, "error_formula")
  vars <- unique(c(all.vars(stats::delete.response(stats::terms(carrier))), all.vars(ef)))
  if (!all(vars %in% names(mf)) || !all(vapply(mf[vars], is.factor, logical(1)))) return(NULL)
  key <- function(d) do.call(paste, c(lapply(d[vars], as.character), sep = "\r"))
  # (a) rows of the data with the response missing
  cl <- attr(aovl, "call")
  dat <- tryCatch(eval(cl$data, environment(stats::formula(stats::terms(aovl)))), error = function(e) NULL)
  yv <- all.vars(stats::formula(stats::terms(aovl))[[2]])
  # (rows a 'subset' left out are not missing values: then only (b))
  if (is.null(cl$subset) && is.data.frame(dat) && all(c(vars, yv) %in% names(dat))) {
    keep <- stats::complete.cases(dat[vars])
    miss <- keep & !stats::complete.cases(dat[yv])
    if (any(miss)) {
      full <- dat[keep, vars, drop = FALSE]
      full[] <- lapply(names(full), function(v) factor(full[[v]], levels = levels(mf[[v]])))
      if (anyNA(full) || anyDuplicated(key(full))) return(NULL)
      return(list(full = full, obs = match(key(full), key(mf))))
    }
  }
  # (b) the crossing of the design factors, when only a few of its cells are missing
  full <- do.call(expand.grid, c(lapply(mf[vars], function(v) levels(droplevels(v))), KEEP.OUT.ATTRS = FALSE))
  nmiss <- nrow(full) - nrow(mf)
  if (nmiss <= 0 || nmiss > 0.1 * nrow(full) || anyDuplicated(key(mf)) || !all(key(mf) %in% key(full))) return(NULL)
  list(full = full, obs = match(key(full), key(mf)))
}

## The strata of a design (rows of Q' by stratum), as aov() builds them.
pm_aov_qr <- function(ef, dat) {
  op <- options(contrasts = c("contr.helmert", "contr.poly"))
  on.exit(options(op))
  Xe <- stats::model.matrix(stats::as.formula(paste("~", paste(deparse(ef), collapse = ""))), dat)
  qe <- qr(Xe)
  rk <- qe$rank
  asg <- attr(Xe, "assign")[qe$pivot[seq_len(rk)]]
  list(qe = qe, grp = c(asg, rep(max(asg) + 1L, nrow(Xe) - rk)))
}

## Each stratum of a design (pm_aov_qr()): its rows of Q', an orthonormal
## basis U of the treatment design X in them, and its residual df; the
## residual projector of the stratum is I - U U' on its rows.
pm_aov_bases <- function(qe, grp, X) {
  tX <- qr.qty(qe, X)
  lapply(sort(unique(grp)), function(g) {
    sel <- which(grp == g)
    Xi <- tX[sel, , drop = FALSE]
    cols <- colSums(Xi^2) > 1e-8 * max(1, colSums(X^2))
    U <- if (any(cols)) {
      qi <- qr(Xi[, cols, drop = FALSE])
      qr.Q(qi)[, seq_len(qi$rank), drop = FALSE]
    } else matrix(0, length(sel), 0)
    list(sel = sel, U = U, df = length(sel) - ncol(U))
  })
}

## TRUE when every treatment term lies in a single stratum (an orthogonal design).
pm_aov_orthogonal <- function(qe, grp, carrier, dat) {
  op <- options(contrasts = c("contr.helmert", "contr.poly"))
  on.exit(options(op))
  Xh <- stats::model.matrix(stats::delete.response(stats::terms(carrier)), dat)
  qh <- qr(Xh)
  if (qh$rank < ncol(Xh)) return(FALSE)
  Qo <- qr.Q(qh)
  tq <- qr.qty(qe, Qo)
  asg <- attr(Xh, "assign")
  for (k in unique(asg)) {
    e <- rowsum(tq[, asg == k, drop = FALSE]^2, grp)          # strata x columns
    if (any(apply(e, 2, function(v) max(v) < (1 - 1e-8) * sum(v)))) return(FALSE)
    if (length(unique(apply(e, 2, which.max))) > 1) return(FALSE)
  }
  TRUE
}

## The analysis with the missing values estimated (NULL when the completed
## design is not orthogonal). Q (the strata's rotation, from the QR of the
## Error() model) is applied to the few columns needed, never formed: the cost
## grows with n x (rank of the Error() model) x (coefficients + missing
## values), not n^3.
pm_aov_missing <- function(carrier, comp, X, y) {
  if (ncol(X) != length(stats::coef(carrier))) return(NULL)      # aliased in the observed data
  full <- comp$full
  # rows of the carrier's frame (the fixed part) for the observed units
  was_obs <- !is.na(comp$obs)
  rn_aov <- rownames(stats::model.frame(attr(carrier, "aovlist")))
  comp$obs <- match(rn_aov[comp$obs], rownames(stats::model.frame(carrier)))
  if (anyNA(comp$obs[was_obs])) return(NULL)                     # an observed unit not in the fixed part
  ef <- attr(carrier, "error_formula")
  sq <- pm_aov_qr(ef, full)
  if (!pm_aov_orthogonal(sq$qe, sq$grp, carrier, full)) return(NULL)
  qe <- sq$qe
  n <- nrow(full)
  obs <- !is.na(comp$obs)
  m <- which(!obs)
  o <- which(obs)
  no <- length(o)
  nm_ <- length(m)
  Xf <- stats::model.matrix(stats::delete.response(stats::terms(carrier)), full,
                            contrasts.arg = carrier$contrasts)[, colnames(X), drop = FALSE]
  p <- ncol(Xf)
  grp <- sq$grp
  st <- sort(unique(grp))
  tX <- qr.qty(qe, Xf)                                           # Q'X, by stratum rows
  S <- pm_aov_bases(qe, grp, Xf)
  res <- function(i, Z) Z - S[[i]]$U %*% crossprod(S[[i]]$U, Z)
  # rows m of the residual projector of each stratum, R_i[m, ] = Q_i (I - U U') Q_i'[, m]
  Em <- matrix(0, n, nm_)
  Em[cbind(m, seq_len(nm_))] <- 1
  Tm <- qr.qty(qe, Em)                                           # Q'[, m]
  Rm <- lapply(seq_along(st), function(i) {
    Z <- matrix(0, n, nm_)
    Z[S[[i]]$sel, ] <- res(i, Tm[S[[i]]$sel, , drop = FALSE])
    t(qr.qy(qe, Z))                                               # nm_ x n
  })
  # the missing values, stratum by stratum from the lowest: y_m = C y_o + N z
  C <- matrix(0, nm_, no)
  N <- diag(nm_)
  red <- integer(length(st))
  for (k in rev(seq_along(st))) {
    if (!ncol(N) || S[[k]]$df <= 0) next
    Amm <- Rm[[k]][, m, drop = FALSE]
    Amo <- Rm[[k]][, o, drop = FALSE]
    SS <- crossprod(N, Amm %*% N)
    es <- eigen((SS + t(SS)) / 2, symmetric = TRUE)
    pos <- es$values > 1e-8 * max(1, abs(es$values))
    if (!any(pos)) next
    Sinv <- es$vectors[, pos, drop = FALSE] %*% (t(es$vectors[, pos, drop = FALSE]) / es$values[pos])
    C <- C - N %*% Sinv %*% crossprod(N, Amo + Amm %*% C)
    red[k] <- sum(pos)
    N <- N %*% es$vectors[, !pos, drop = FALSE]
  }
  yo <- y[comp$obs[o]]
  yc <- numeric(n)
  yc[o] <- yo
  yc[m] <- drop(C %*% yo)
  # stratum mean squares of the completed data, residual df less the values each determined
  qyc <- qr.qty(qe, yc)
  dfs <- vapply(S, `[[`, numeric(1), "df") - red
  ss <- vapply(seq_along(st), function(i) sum(res(i, qyc[S[[i]]$sel])^2), numeric(1))
  ms <- ifelse(dfs > 0, ss / pmax(dfs, 1), NA_real_)
  known <- which(dfs > 0)
  if (!length(known)) return(NULL)
  donor <- vapply(seq_along(st), function(i) {
    if (dfs[i] > 0) return(i)
    below <- known[known > i]
    if (length(below)) min(below) else max(known[known < i])
  }, integer(1))
  xi <- ms[donor]
  # coefficients of the completed design (each term from its own stratum),
  # b = Bc y_c with Bc = W^-1 sum_i X_i'Q_i'/xi_i, and in the observed data
  # B = Bc E (E: y_o -> y_c)
  W <- Reduce(`+`, lapply(seq_along(st), function(i) crossprod(tX[S[[i]]$sel, , drop = FALSE]) / xi[i]))
  Ztx <- tX
  for (i in seq_along(st)) Ztx[S[[i]]$sel, ] <- tX[S[[i]]$sel, , drop = FALSE] / xi[i]
  Bc <- solve(W, t(qr.qy(qe, Ztx)))                              # p x n
  B <- Bc[, o, drop = FALSE] + Bc[, m, drop = FALSE] %*% C        # p x no
  b <- drop(B %*% yo)
  # exact covariance: sum_d ms_d B Q[o, strata of d] Q[o, strata of d]' B'
  Zb <- matrix(0, n, p)
  Zb[o, ] <- t(B)
  QZ <- qr.qty(qe, Zb)                                            # Q[o, ]' B'
  Jac <- lapply(known, function(d) {
    rows <- unlist(lapply(which(donor == d), function(i) S[[i]]$sel))
    crossprod(QZ[rows, , drop = FALSE])
  })
  V <- Reduce(`+`, Map(function(J, d) J * ms[d], Jac, known))
  dimnames(V) <- list(colnames(X), colnames(X))
  names(b) <- colnames(X)
  parts <- structure(list(vcov_beta = V, vcov_varpar = diag(2 * ms[known]^2 / dfs[known], nrow = length(known)),
                          Jac_list = Jac), class = "pm_satt")
  nms <- c("(Intercept)", attr(stats::terms(stats::as.formula(paste("~", paste(deparse(ef), collapse = "")))),
                               "term.labels"), "Within")[st + 1L]
  message("aov() with Error() and ", nm_, " missing value", if (nm_ > 1) "s", ": estimated as in Genstat's ANOVA ",
          "(the residual sum of squares of the lowest stratum made smallest) and the design analysed by strata, ",
          "with exact standard errors. REML (e.g. lmer() with the Error() terms as random effects) uses the ",
          "information in all strata and is the better analysis when plots are missing.")
  mp <- structure(list(coef = b, vcov = V, df = NA_real_, estimable = rep(TRUE, length(b))), class = "mymodelparm")
  list(mp = mp, parts = parts,
       strata = structure(data.frame(stratum = nms, df = dfs, ms = ms),
                          missing = data.frame(full[m, , drop = FALSE], estimate = yc[m], row.names = NULL)))
}
