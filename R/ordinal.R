## Cumulative link models: MASS::polr(), ordinal::clm() and ordinal::clmm()
## (decisions D1-D6 of dev/review-classes.md, 2026-11-25). clmm: the fixed part,
## random effects at 0 (a typical subject, as emmeans; D4); its covariance of
## (b, alpha) is the block of vcov() without the random-effect parameters.

## The classes of a class response (one list for every check in the package);
## multinom and the vglm fits of pm_cat_is() go to R/categorical.R.
pm_ord_classes <- c("polr", "clm", "clmm", "vglm", "multinom")
##
## Both fit P(Y <= j) = F(theta_j - eta), eta = x'b without an intercept, cut
## points theta = J alpha (J = I for polr and for clm's flexible thresholds;
## clm's tJac otherwise). The means are on the latent scale, x'b - mean(theta)
## (as emmeans' mode = "latent"): the coefficients are given as
## ((Intercept) = -mean(theta), b), so the model's own design (with its
## intercept column) gives the means, and the covariance is that of this
## linear map of (b, alpha). A difference of latent means (logit link) is the
## log odds ratio of a higher class, log[P(Y > j) / P(Y <= j)], the same for
## every cut point j. Wald tests are z tests (pm_df_rule(), as fixed-dispersion
## GLMs). Aliased coefficients are NA in pm_coef() and left out of pm_vcov()
## (the adapter contract, ?pm_coef).
##
## VGAM::vglm (B3a, dev/review-vglm.md): ordinal families with parallel = TRUE
## (cumulative, propodds, acat, cratio, sratio). Their predictors are
## eta_j = alpha_j + x'b; the sign of b differs between families and with
## 'reverse', so it is found from the fit (s: +1 when a larger eta makes higher
## classes more likely), and the latent mean is s (x'b + mean(alpha)), equal to
## clm's x'b - mean(theta) for the same cumulative model. Class probabilities
## come from the family's own linkinv, their gradients by finite differences.

## b (aliased coefficients NA), alpha, J, the joint covariance of the
## estimable b and alpha (in that order), the link and the class labels of a
## polr or clm fit; stops for what is not handled.
pm_ord_parts <- function(model) {
  if (pm_cat_is(model)) pm_cat_stop(model)
  if (inherits(model, "vglm")) return(pm_ord_parts_vglm(model))
  if (inherits(model, "polr")) {
    if (is.null(model$Hessian)) {
      stop("The polr fit has no Hessian: refit with polr(..., Hess = TRUE).", call. = FALSE)
    }
    b <- stats::coef(model)
    alpha <- model$zeta
    J <- diag(length(alpha))
    link <- model$method
    lev <- model$lev
  } else {
    if (inherits(model, "clmm") && !is.null(model$call$nominal)) {
      stop("clmm() does not fit nominal effects: it ignored 'nominal ='; refit without it.", call. = FALSE)
    }
    pm_clm_link_check(model)
    b <- model$beta
    alpha <- model$alpha
    J <- model$tJac
    link <- model$link
    lev <- model$y.levels
  }
  V <- pm_ord_vcov(model)
  nm <- c(names(b)[!is.na(b)], names(alpha))
  if (is.null(V) || anyNA(V) || !all(nm %in% rownames(V))) {
    stop("The ", class(model)[1], " fit has no usable covariance matrix of its coefficients.", call. = FALSE)
  }
  list(b = b, alpha = alpha, J = J, V = V[nm, nm, drop = FALSE], link = link, levels = lev, sb = 1, st = -1,
       family = "cumulative", probs = NULL)
}

## vcov() of an ordinal fit. A clmm whose random-effect SD is estimated at 0
## on ordinal's log-SD route has a Hessian with a zero row (the SD), and
## vcov() stops; the covariance of the other parameters is then the inverse
## of their block, the SD treated as known (as ordinal does on its other
## route, where SDs at 0 are left out); attribute "sd_fixed".
pm_ord_vcov <- function(model) {
  V <- tryCatch(stats::vcov(model), error = function(e) e)
  if (!inherits(V, "error")) return(V)
  sd <- if (inherits(model, "clmm")) vapply(model$ST, function(m) as.matrix(m)[1, 1], numeric(1))
  H <- model$Hessian
  st <- grepl("^ST[0-9]+$", rownames(H))
  if (is.null(sd) || !any(sd <= 1e-3) || !any(st)) stop(conditionMessage(V), call. = FALSE)
  out <- tryCatch(solve(H[!st, !st, drop = FALSE]), error = function(e) NULL)
  if (is.null(out)) stop(conditionMessage(V), call. = FALSE)
  structure(out, sd_fixed = TRUE)
}

## Links with an estimated parameter stop (clm).
pm_clm_link_check <- function(model) {
  if (!is.null(model$lambda) || model$link %in% c("Aranda-Ordaz", "log-gamma")) {
    stop("The ", model$link, " link has an estimated parameter (lambda), which the standard errors would treat as ",
         "known; use a logit, probit, cloglog, loglog or cauchit link.", call. = FALSE)
  }
}

## The parts of a vglm fit (B3a): ordinal families with parallel = TRUE only.
## eta = J alpha + x'b (J: the constraint matrix of the intercepts, I for free
## cut points, e.g. cbind(1, 0:3) for equidistant ones). probs(b, alpha, x):
## the family's linkinv; deta(eta): its Jacobian in eta (finite differences
## on eta, with steps scaled to eta, then the chain rule in pm_ord_class()).
pm_ord_parts_vglm <- function(model) {
  fi <- pm_ord_vglm_info(model)
  cm <- VGAM::constraints(model)
  pm_vglm_smart_stop(model)
  J <- cm[["(Intercept)"]]
  cf <- VGAM::coef(model)
  ic <- grep("^\\(Intercept\\)", names(cf))
  alpha <- cf[ic]
  b <- cf[-ic]
  V <- VGAM::vcov(model)
  nm <- c(names(b), names(alpha))
  if (is.null(V) || anyNA(V) || !all(nm %in% rownames(V)) || nrow(J) < 1 || ncol(J) != length(alpha)) {
    stop("The vglm fit has no usable covariance matrix of its coefficients.", call. = FALSE)
  }
  # (two rows: VGAM's acat(reverse = TRUE) linkinv drops a one-row matrix)
  linkinv <- function(eta) as.vector(model@family@linkinv(rbind(eta, eta), extra = model@extra)[1, ])
  lev <- fi$levels
  # the sign: does a larger eta (all predictors) raise the mean class? At the
  # average linear predictor of the data (the cut points alone, x = 0, can be
  # far from the data, where the probabilities are flat)
  mc <- function(e) sum(seq_along(lev) * linkinv(e))
  e0 <- colMeans(model@predictors)
  s <- sign(mc(e0 + 0.01) - mc(e0 - 0.01))
  if (s == 0) stop("Could not find the direction of the vglm family's predictors.", call. = FALSE)
  list(b = b, alpha = alpha, J = J, V = V[nm, nm, drop = FALSE], link = pm_ord_vglm_link(fi), levels = lev,
       sb = s, st = s, family = fi$family,
       probs = function(b, alpha, x) linkinv(as.vector(J %*% alpha) + sum(x * b)),
       deta = function(eta) pm_fd_jac(linkinv, eta))
}

## Family, link and class labels of a vglm fit, without the covariance (for
## labels); stops for families that are not ordinal.
pm_ord_vglm_info <- function(model) {
  fam <- model@family@vfamily[1]
  if (pm_cat_is(model)) pm_cat_stop(model)
  if (!fam %in% c("cumulative", "acat", "cratio", "sratio")) {
    stop("predictmeans handles vglm fits of the families cumulative, propodds, acat, cratio, sratio and ",
         "multinomial; not ", fam, " (vglm fits of other families are not supported).", call. = FALSE)
  }
  list(family = fam, link = unname(sub("link$", "", model@misc$link[1])), reverse = isTRUE(model@misc$reverse),
       levels = colnames(model@y))
}

## The link in clm's convention (P(Y <= j) = F(theta_j - eta)): a cumulative
## fit with reverse = TRUE models P(Y > j), so its cloglog is clm's loglog and
## the other way round (symmetric links are the same).
pm_ord_vglm_link <- function(fi) {
  if (fi$family == "cumulative" && fi$reverse) {
    return(switch(fi$link, cloglog = "loglog", loglog = "cloglog", fi$link))
  }
  fi$link
}

## The checks of an ordinal fit and of the options of a call, for every front
## end (predictmeans(), covariatemeans()): pm_ord_parts() stops for fits not
## handled; latent means have no back-transformation and no residual
## variance, so 'trans', 'bias_adjust', 'sigma' (and 'vcov', 'permlist', not
## built for the latent parametrisation) stop: the class scale is mean_type's.
pm_ord_check <- function(model, trans = NULL, bias_adjust = FALSE, sigma = NULL, vcov = NULL, permlist = NULL) {
  pm_ord_parts(model)
  pm_class_opt_stop(model, paste0("the means are on the latent scale; for the class scale use mean_type = \"prob\" ",
                                  "or \"mean_class\" in predictmeans()",
                                  if (inherits(model, "clmm") && !isFALSE(bias_adjust))
                                    ", where bias_adjust = TRUE gives population-average class probabilities"),
                    trans, bias_adjust, sigma, vcov, permlist)
  invisible(TRUE)
}

## The options built for one scale with a residual variance, which models of
## a class response do not take: a stop naming those given ('why': the scale
## of the means). One place for pm_ord_check() and pm_cat_check().
pm_class_opt_stop <- function(model, why, trans = NULL, bias_adjust = FALSE, sigma = NULL, vcov = NULL, permlist = NULL) {
  bad <- c(trans = !is.null(trans) && !identical(trans, I), bias_adjust = !isFALSE(bias_adjust),
           sigma = !is.null(sigma), vcov = !is.null(vcov), permlist = !is.null(permlist))
  if (any(bad)) {
    stop("'", paste(names(bad)[bad], collapse = "', '"), "' cannot be used with ", class(model)[1], " models: ",
         why, ".", call. = FALSE)
  }
}

## VGAM's smart terms are not rebuilt on the grid of the means: a stop.
pm_vglm_smart_stop <- function(model) {
  if (any(grepl("^sm\\.", all.names(stats::formula(model))))) {
    stop("VGAM's smart terms (sm.bs, sm.ns, sm.poly, sm.scale) cannot be rebuilt on the grid of the means; refit ",
         "with bs(), ns(), poly() or scale().", call. = FALSE)
  }
}

## The coefficients on the latent scale and their covariance: ((Intercept) =
## st mean(theta), sb b) (polr, clm: -mean(theta), b; vglm: s (mean(alpha), b)),
## A (b, alpha) with A the linear map (estimable b only).
pm_ord_latent <- function(model) {
  p <- pm_ord_parts(model)
  be <- p$b[!is.na(p$b)]
  nb <- length(be)
  A <- rbind(c(rep(0, nb), p$st * colMeans(p$J)), cbind(p$sb * diag(nb), matrix(0, nb, length(p$alpha))))
  icpt <- sum(A[1, ] * c(be, p$alpha))
  V <- A %*% p$V %*% t(A)
  dimnames(V) <- list(c("(Intercept)", names(be)), c("(Intercept)", names(be)))
  list(coef = c("(Intercept)" = icpt, p$sb * p$b), vcov = V)
}

#' @exportS3Method pm_coef polr
pm_coef.polr <- function(model, ...) pm_ord_latent(model)$coef
#' @exportS3Method pm_coef clm
pm_coef.clm <- function(model, ...) pm_ord_latent(model)$coef
#' @exportS3Method pm_vcov polr
pm_vcov.polr <- function(model, ...) pm_ord_latent(model)$vcov
#' @exportS3Method pm_vcov clm
pm_vcov.clm <- function(model, ...) pm_ord_latent(model)$vcov
#' @exportS3Method pm_coef clmm
pm_coef.clmm <- pm_coef.clm
#' @exportS3Method pm_vcov clmm
pm_vcov.clmm <- pm_vcov.clm
## no back-transformation of the latent means (class probabilities: mean_type = "prob")
## the contrasts the fit used (clmm has no model.matrix() method: the default
## would rebuild the design with the current options)
#' @exportS3Method pm_contrasts clmm
pm_contrasts.clmm <- function(model, ...) model$contrasts
#' @exportS3Method pm_linkinv polr
pm_linkinv.polr <- function(model, ...) NULL
#' @exportS3Method pm_linkinv clm
pm_linkinv.clm <- function(model, ...) NULL
#' @exportS3Method pm_linkinv clmm
pm_linkinv.clmm <- function(model, ...) NULL
#' @exportS3Method pm_coef vglm
pm_coef.vglm <- function(model, ...) pm_ord_latent(model)$coef
#' @exportS3Method pm_vcov vglm
pm_vcov.vglm <- function(model, ...) pm_ord_latent(model)$vcov
#' @exportS3Method pm_linkinv vglm
pm_linkinv.vglm <- function(model, ...) NULL
#' @exportS3Method pm_contrasts vglm
pm_contrasts.vglm <- function(model, ...) model@contrasts

## z tests: the latent scale has no dispersion to estimate (D3; as emmeans and
## the fixed-dispersion GLM rule, dev/sim-df-glm/).
#' @exportS3Method NULL
pm_df_rule.polr <- function(model, df_method, kr_info) {
  if (df_method != "auto") {
    warning("'df_method = \"", df_method, "\"' does not apply to ", class(model)[1], " models; ignored.", call. = FALSE)
  }
  list(df_type = "asymptotic", vcov_adj = NULL)
}
#' @exportS3Method NULL
pm_df_rule.clm <- pm_df_rule.polr
## clmm: inner-outer df by default (pm_df_method(), as GLMMs: z rejected 11% of
## true nulls at 5% with 6 clusters for a between-cluster term); this rule
## (z) is the fallback when they cannot be found.
#' @exportS3Method NULL
pm_df_rule.clmm <- pm_df_rule.polr
#' @exportS3Method NULL
pm_df_rule.vglm <- pm_df_rule.polr

## Class probabilities P(Y = j) (mean_type = "prob") or the mean class
## sum_j j P(Y = j) (mean_type = "mean_class") at each mean: rows of K (the
## latent means; averages over other factors are taken on the latent scale,
## D2), with delta-method SEs from the covariance of (b, alpha). Limits:
## probabilities on the logit scale (they stay in (0, 1)), the mean class on
## its own scale. emmeans' mode = "prob" / "mean.class" give the same values
## for the cells of the model (they average probabilities over other factors).
## The distribution function F of a cumulative link and its density.
pm_ord_link_fd <- function(link) {
  switch(link, logistic = , logit = list(stats::plogis, stats::dlogis),
         probit = list(stats::pnorm, stats::dnorm),
         cloglog = list(function(x) 1 - exp(-exp(x)), function(x) exp(x - exp(x))),
         loglog = list(function(x) exp(-exp(-x)), function(x) exp(-x - exp(-x))),
         cauchit = list(stats::pcauchy, stats::dcauchy),
         stop("link '", link, "' is not supported.", call. = FALSE))
}

## Population-average class probabilities of a clmm (B4 G2, dev/review-b4.md):
## P(Y <= j | x) = E F(z - s U), z = theta_j - x'b, U ~ N(0, 1), s^2 the sum of
## the random-intercept variances (z' Sigma z for intercepts: one normal
## number, so one integral; pm_ord_pop_int(): a 60-node Gauss-Hermite rule for
## s <= 0.5, integrate() split where z - s u = 0 above, where a fixed rule was
## off by up to 5e-3 at s = 5 for cloglog), with its derivatives in z and s. Returns s, its gradient in the
## random-effect parameters of vcov() (named), which are log SDs for one
## random intercept fitted the default way and SDs otherwise (decided from the
## optimiser's values, model$optRes$par); SDs at 0 (below 1e-3, left out of
## vcov() by ordinal) are fixed; V, the covariance of (alpha, b, those
## parameters). Random slopes (any term other than an intercept) stop: the
## variance then depends on the covariate values.
pm_ord_pop <- function(model) {
  if (!pm_ord_intercepts_only(model)) {
    stop("Population-average probabilities (bias_adjust = TRUE) are available for random intercepts; this clmm ",
         "has random slopes, whose variance depends on the covariate values. Use bias_adjust = FALSE (a typical ",
         "subject).", call. = FALSE)
  }
  sd <- vapply(model$ST, function(m) as.matrix(m)[1, 1], numeric(1))
  s <- sqrt(sum(sd^2))
  par <- utils::tail(as.vector(model$optRes$par), length(sd))
  logsd <- length(par) == length(sd) && all(abs(par - log(pmax(sd, 1e-300))) < abs(par - sd))
  V <- pm_ord_vcov(model)
  free <- sd > 1e-3 | (logsd & !isTRUE(attr(V, "sd_fixed")))
  stn <- grep("^ST[0-9]+$", rownames(V), value = TRUE)
  if (length(stn) != sum(free)) stop("internal: the random-effect parameters of the clmm fit are not matched.",
                                     call. = FALSE)
  ds <- if (s > 0) (sd / s * if (logsd) sd else 1)[free] else rep(0, sum(free))
  list(s = s, ds = stats::setNames(ds, stn), V = V)
}

## TRUE when every random term of a clmm is an intercept, (1 | g) (a scalar
## term (0 + x | g) has a 1 x 1 ST too, but its variance depends on x).
pm_ord_intercepts_only <- function(model) {
  f <- stats::formula(model)
  bars <- pm_plain_bars(f, pm_special_terms(f))
  all(vapply(bars, function(b) identical(deparse(b[[2]]), "1"), logical(1)))
}

## Gauss-Hermite rule for a standard normal (Golub-Welsch): E g(U) = sum w g(x).
pm_gh_rule <- function(n = 60) {
  k <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1)] <- J[cbind(k + 1, k)] <- sqrt(k)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}

## E F(z - s U), E f(z - s U) and d/ds E F(z - s U) = -E[U f(z - s U)] for
## each z (Fd: the link's distribution function and density).
pm_ord_pop_int <- function(z, s, Fd) {
  if (s <= 0) return(list(F = Fd[[1]](z), f = Fd[[2]](z), ds = rep(0, length(z))))
  if (s <= 0.5) {
    # a small SD: the integrand is smooth, and a 60-node Gauss-Hermite rule is accurate to 2e-13 for every link
    # (integrate() can stop with "roundoff error" there)
    gh <- pm_gh_rule(60)
    A <- outer(z, s * gh$x, "-")
    return(list(F = drop(Fd[[1]](A) %*% gh$w), f = drop(Fd[[2]](A) %*% gh$w),
                ds = -drop(Fd[[2]](A) %*% (gh$w * gh$x))))
  }
  # over |u| <= 12 (the normal tail beyond is below 1e-32), split where z - s u = 0 when that is inside (the
  # step of F and the peak of f; with a small s it is outside and the integrand is smooth)
  one <- function(g, z0) {
    u0 <- z0 / s
    br <- c(-12, if (abs(u0) < 12) u0, 12)
    sum(vapply(seq_len(length(br) - 1), function(i) stats::integrate(function(u) g(u) * stats::dnorm(u), br[i], br[i + 1],
                                                                     rel.tol = 1e-11, abs.tol = 1e-15,
                                                                     subdivisions = 500L)$value, 0))
  }
  out <- lapply(z, function(z0) c(one(function(u) Fd[[1]](z0 - s * u), z0), one(function(u) Fd[[2]](z0 - s * u), z0),
                                   -one(function(u) u * Fd[[2]](z0 - s * u), z0)))
  out <- matrix(unlist(out), 3)
  list(F = out[1, ], f = out[2, ], ds = out[3, ])
}

pm_ord_class <- function(model, K, what = c("prob", "mean_class"), conf_level = 0.95, pop = FALSE) {
  what <- match.arg(what)
  p <- pm_ord_parts(model)
  if (pop) {                                     # population average: (b, alpha) and the random-effect SDs
    pp <- pm_ord_pop(model)
    nm <- c(rownames(p$V), names(pp$ds))
    p$V <- pp$V[nm, nm, drop = FALSE]
  }
  Fd <- if (is.null(p$probs)) pm_ord_link_fd(p$link)
  be <- p$b[!is.na(p$b)]
  X <- as.matrix(K)[, names(be), drop = FALSE]
  th <- as.vector(p$J %*% p$alpha)
  nj <- length(p$levels) - 1
  z <- pm_ci_crit(1 - conf_level, Inf)                          # z: the df of these models (pm_df_rule)
  rows <- lapply(seq_len(nrow(X)), function(k) {
    x <- X[k, ]
    if (anyNA(x)) return(NULL)
    if (is.null(p$probs)) {
      zz <- th - sum(x * be)
      Gs <- NULL
      if (pop) {                                 # E F(zz - s U), its derivatives in zz and in s
        I <- pm_ord_pop_int(zz, pp$s, Fd)
        f <- I$f
        cum <- c(0, I$F, 1)
        Gs <- outer(I$ds, pp$ds)
      } else {
        f <- Fd[[2]](zz)
        cum <- c(0, Fd[[1]](zz), 1)
      }
      # gradients of P(Y <= j) with respect to (b, alpha[, random-effect SDs]), j = 0, ..., J
      G <- rbind(0, cbind(-outer(f, x), f * p$J, Gs), 0)
      P <- diff(cum)
      GP <- G[-1, , drop = FALSE] - G[-(nj + 2), , drop = FALSE]
    } else {                                     # vglm: the family's linkinv; d eta / d(b, alpha) = [x 1', J]
      eta <- as.vector(p$J %*% p$alpha) + sum(x * be)
      P <- p$probs(be, p$alpha, x)
      De <- p$deta(eta)
      GP <- De %*% cbind(matrix(x, length(eta), length(x), byrow = TRUE), p$J)
    }
    if (what == "mean_class") {
      g <- colSums(seq_along(P) * GP)
      m <- sum(seq_along(P) * P)
      se <- sqrt(sum((g %*% p$V) * g))
      # limits on the logit of (m - 1) / (J - 1), so they stay within the classes 1 to J
      u <- (m - 1) / nj
      lu <- stats::qlogis(u)
      su <- se / nj / (u * (1 - u))
      return(data.frame(row = k, estimate = m, std.error = se, conf.low = 1 + nj * stats::plogis(lu - z * su),
                        conf.high = 1 + nj * stats::plogis(lu + z * su)))
    }
    se <- sqrt(rowSums((GP %*% p$V) * GP))
    lg <- stats::qlogis(P)
    sl <- se / (P * (1 - P))                       # SE of logit(P)
    data.frame(row = k, class = factor(p$levels, levels = p$levels), estimate = P, std.error = se,
               conf.low = stats::plogis(lg - z * sl), conf.high = stats::plogis(lg + z * sl))
  })
  structure(do.call(rbind, rows), s = if (pop) pp$s)
}

## Add the class table to a predictmeans() result: $probs, rows in the order
## of the means (matched by their labels, not position: atvar re-sorts them),
## with the factor columns of modelterm in front.
pm_ord_attach <- function(res, model, K, what, conf_level, pop = FALSE) {
  means <- res$means
  vars <- intersect(strsplit(res$settings$modelterm, ":")[[1]], names(means))
  key <- do.call(paste, c(lapply(means[vars], as.character), sep = ":"))
  rows <- match(key, rownames(K))
  if (anyNA(rows)) stop("internal: the means do not match the rows of K.", call. = FALSE)
  tab <- pm_ord_class(model, K[rows, , drop = FALSE], what, conf_level, pop)
  out <- cbind(means[tab$row, vars, drop = FALSE], tab[setdiff(names(tab), "row")])
  rownames(out) <- NULL
  attr(out, "what") <- what
  attr(out, "population_average") <- pop
  res$probs <- out
  if (pop) res$settings$bias_adjust <- attr(tab, "s")           # (the SD averaged over, as for other models)
  pm_add_advice(res, pm_ordinal_advice(model, res$settings$modelterm, what, pop))
}

## What a difference of latent means is, for the summary line: the scale of
## the family and link (a positive difference: higher classes more likely).
pm_ord_scale_text <- function(family, link) {
  fam <- if (is.null(family)) "cumulative" else family
  if (identical(link, "logit") && fam == "cumulative") return("log odds ratios of a higher class, the same at every cut point")
  what <- switch(fam, acat = "adjacent-category", cratio = , sratio = "continuation-ratio", "cumulative")
  paste0(what, " ", if (link %in% c("logit", "log")) "log odds ratios" else paste0("differences on the ", link, " scale"),
         ", the same at every cut point (a positive difference: higher classes more likely)")
}
