## Categorical responses without a single latent scale (B3b, dev/review-vglm.md
## E2): nnet::multinom(), VGAM::vglm() with multinomial(), and vglm ordinal
## families (cumulative, acat, cratio, sratio) whose slopes are not parallel,
## and ordinal::clm() with nominal or scale effects (B4 G1, dev/review-b4.md).
## Each of the M = J - 1 linear predictors has its own coefficients, so a
## mean is a vector. predictmeans() gives, for each mean of 'modelterm':
##   mean_type = "prob" (the default)  the probability of each class: averaged
##        over the other factors on the probability scale (as emmeans'
##        mode = "prob"), compared within each class (differences of
##        probabilities, delta-method SEs, z tests);
##   mean_type = "link"   each linear predictor (the model's own scale: exact
##        Wald tests; for multinomial models each depends on the reference
##        class);
##   mean_type = "mean_class"   sum_j j P(Y = j), ordinal families only.
## The class (or predictor) is a factor of the means, and the pairs are
## within its levels (an 'atvar'), so the pairs, letters, LSDs and plots of
## predictmeans() work unchanged on the stacked means, as for the overall
## means of zero-inflated models (pm_zi_overall_est()): K the identity, the
## coefficients the means, their covariance D V D' (D: the gradient in the
## coefficients theta of the fit).
##
## Kmatrix() builds the grid and the averaging weights from the design of
## the fit (one column per term column, x), through a carrier (class
## "pm_catfit", pm_cat_carrier()) made by predictmeans() only; the fit itself
## stops in the other functions (pm_coef()), which do not handle a vector mean.

## TRUE for the fits handled here.
pm_cat_is <- function(model) {
  if (inherits(model, "multinom")) return(TRUE)
  if (inherits(model, "clm")) return(length(model$nom.terms) > 0 || length(model$S.terms) > 0 || length(model$zeta) > 0)
  if (!inherits(model, "vglm")) return(FALSE)
  fam <- model@family@vfamily[1]
  fam == "multinomial" || (fam %in% c("cumulative", "acat", "cratio", "sratio") && !pm_ord_vglm_parallel(model))
}

## TRUE when every slope of a vglm has one coefficient for all predictors
## (parallel = TRUE; its constraint matrix a column of ones).
pm_ord_vglm_parallel <- function(model) {
  cm <- VGAM::constraints(model)
  slope <- cm[names(cm) != "(Intercept)"]
  all(vapply(slope, function(C) ncol(C) == 1 && all(C == 1), logical(1)))
}

## The parts of a fit: theta (the coefficients, named as vcov()), V their
## covariance, the design columns 'xnames', and the map to the predictors:
## eta(x) = H (theta * x[col]) (H: M x q; col: the design column of each
## coefficient). pd(eta), eta a matrix (rows x M): list(P, the class
## probabilities (rows x J), dP their derivatives in eta (rows x J x M)). levels: the classes; predictors: the names of
## the linear predictors; reference: the reference class of a multinomial
## fit; ordinal: TRUE for the ordinal families. Fits whose predictors are
## not linear in theta (clm scale effects) add 'scale' (pm_clm_cat_parts()),
## and fits with several design matrices 'basis' (the design columns of
## newdata).
pm_cat_parts <- function(model) {
  if (inherits(model, "clm")) return(pm_clm_cat_parts(model))
  if (inherits(model, "multinom")) {
    if (isTRUE(model$censored) || !isTRUE(model$softmax) && !isTRUE(model$entropy)) {
      stop("predictmeans handles multinom fits of a class response (softmax, or a two-class response); not censored ",
           "or other nnet fits.", call. = FALSE)
    }
    if (isTRUE(model$decay > 0)) {
      stop("The multinom fit is penalised (decay = ", model$decay, "): its Hessian is not the information of the ",
           "estimates, so SEs and tests would be wrong; refit with decay = 0.", call. = FALSE)
    }
    # nnet keeps arbitrary values for aliased columns (and vcov() a g-inverse): the means would not be unique
    X <- tryCatch(stats::model.matrix(stats::terms(model), stats::model.frame(model), contrasts.arg = model$contrasts),
                  error = function(e) NULL)
    qx <- if (!is.null(X)) qr(X)
    if (!is.null(qx) && qx$rank < ncol(X)) {
      stop("The design of the multinom fit is rank deficient (aliased: ",
           paste(colnames(X)[qx$pivot[-seq_len(qx$rank)]], collapse = ", "), "): nnet keeps arbitrary values for ",
           "these coefficients, so the means would not be unique. Refit with a design of full rank (leave out the ",
           "aliased term, or code nested factors as one factor).", call. = FALSE)
    }
    lev <- if (!is.null(model$lev)) model$lev else model$lab
    B <- stats::coef(model)
    if (!is.matrix(B)) B <- matrix(B, 1, dimnames = list(lev[2], names(B)))
    M <- nrow(B)
    p <- ncol(B)
    theta <- as.vector(t(B))
    names(theta) <- if (M == 1) colnames(B) else paste(rep(rownames(B), each = p), colnames(B), sep = ":")
    V <- tryCatch(stats::vcov(model), error = function(e) NULL)
    H <- kronecker(diag(M), matrix(1, 1, p))
    col <- rep(seq_len(p), M)
    xnames <- colnames(B)
    # softmax, reference class first; d P_i / d eta_m = P_i (delta(i, m + 1) - P_(m + 1))
    pd <- function(eta) {
      E <- cbind(0, eta)
      E <- exp(E - apply(E, 1, max))
      P <- E / rowSums(E)
      dP <- array(0, c(nrow(P), ncol(P), M))
      for (m in seq_len(M)) dP[, , m] <- P * (matrix(seq_len(ncol(P)) == m + 1, nrow(P), ncol(P), byrow = TRUE) - P[, m + 1])
      list(P = P, dP = dP)
    }
    ref <- lev[1]
    predictors <- paste0("log(P[", lev[-1], "]/P[", ref, "])")
    family <- "multinomial"
  } else {
    fam <- model@family@vfamily[1]
    if (!is.null(model@control$xij)) {
      stop("vglm fits with 'xij' (covariates that differ between predictors) are not supported.", call. = FALSE)
    }
    pm_vglm_smart_stop(model)
    cm <- VGAM::constraints(model)
    theta <- VGAM::coef(model)
    V <- tryCatch(VGAM::vcov(model), error = function(e) NULL)
    H <- do.call(cbind, unname(cm))
    col <- rep(seq_along(cm), vapply(cm, ncol, integer(1)))
    xnames <- names(cm)
    if (ncol(H) != length(theta)) stop("internal: the vglm constraints do not match its coefficients.", call. = FALSE)
    # the family's linkinv on all rows at once (two rows at least: some drop a one-row matrix); its
    # derivatives by central differences, one predictor at a time (2M calls)
    li <- function(eta) {
      e <- if (nrow(eta) == 1) rbind(eta, eta) else eta
      as.matrix(model@family@linkinv(e, extra = model@extra))[seq_len(nrow(eta)), , drop = FALSE]
    }
    pd <- function(eta) {
      P <- li(eta)
      dP <- array(0, c(nrow(P), ncol(P), ncol(eta)))
      for (m in seq_len(ncol(eta))) {
        h <- pm_fd_h(eta[, m])
        ep <- em <- eta
        ep[, m] <- eta[, m] + h
        em[, m] <- eta[, m] - h
        dP[, , m] <- (li(ep) - li(em)) / (2 * h)
      }
      list(P = P, dP = dP)
    }
    lev <- colnames(model@y)
    predictors <- model@misc$predictors.names
    ref <- NULL
    if (fam == "multinomial") {                   # (VGAM's "log(mu[,1]/mu[,3])", by level names)
      ref <- lev[model@extra$use.refLevel]
      predictors <- paste0("log(P[", setdiff(lev, ref), "]/P[", ref, "])")
    }
    family <- fam
  }
  if (is.null(V) || anyNA(V) || any(!is.finite(V)) || !all(names(theta) %in% rownames(V))) {
    stop("The ", class(model)[1], " fit has no usable covariance matrix of its coefficients (a rank-deficient design, ",
         "or classes that never occur in some cells?).", call. = FALSE)
  }
  V <- as.matrix(V)[names(theta), names(theta), drop = FALSE]
  list(theta = theta, V = V, H = H, col = col, xnames = xnames, pd = pd, levels = lev,
       predictors = predictors, reference = ref, family = family, ordinal = family != "multinomial")
}

## Means and their gradient in theta. X: rows of the design (columns
## xnames). what = "link": the predictors at the rows of X (linear: X are
## the K rows); "prob" / "mean_class": at the grid rows X, averaged with the
## weights W (rows: the means). Returns list(est, D): est a matrix (means x
## classes, predictors, or one column), D the gradient of as.vector(est)
## (stacked by column: all means of the first class, then the second, ...),
## and, with 'labels' (names of the grid rows), the means that use a grid row
## where the cumulative curves cross ('affected') and where ('crossing').
pm_cat_eval <- function(pt, X, W = NULL, what = c("prob", "link", "mean_class"), labels = NULL) {
  what <- match.arg(what)
  X <- as.matrix(X)[, pt$xnames, drop = FALSE]
  Xc <- X[, pt$col, drop = FALSE]                          # x[col] for each coefficient
  eta <- Xc %*% (pt$theta * t(pt$H))                        # rows x M
  Jm <- lapply(seq_len(nrow(pt$H)), function(m) sweep(Xc, 2, pt$H[m, ], "*"))   # d eta_m / d theta
  if (!is.null(pt$scale)) {
    # clm scale effects: eta / exp(s'zeta); d eta_m / d zeta = -eta_m s
    S <- X[, pt$scale$cols, drop = FALSE]
    e <- exp(-drop(S %*% pt$theta[pt$scale$idx]))
    eta <- eta * e
    Jm <- lapply(seq_along(Jm), function(m) {
      J <- Jm[[m]] * e
      J[, pt$scale$idx] <- -eta[, m] * S
      J
    })
  }
  if (what == "link") return(list(est = eta, D = do.call(rbind, Jm)))
  pd <- pt$pd(eta)
  P <- pd$P
  # non-parallel cumulative fits: the cumulative curves can cross, giving negative "probabilities" at
  # some grid rows; the means that use such a row are not valid (pm_cat_est() makes them NA)
  badg <- if (!is.null(labels)) which(apply(P < -1e-10, 1, any)) else integer(0)
  affected <- if (length(badg)) rowSums(abs(W[, badg, drop = FALSE])) > 0 else rep(FALSE, NROW(W))
  crossing <- if (length(badg)) {
    b <- which(P[badg[1], ] < -1e-10)[1]
    paste0("class '", pt$levels[b], "' has probability ", format(signif(P[badg[1], b], 3)), " at ", labels[badg[1]],
           if (length(badg) > 1) paste0(" (and ", length(badg) - 1, " more grid point", if (length(badg) > 2) "s", ")"))
  }
  # d P_j / d theta = sum_m (d P_j / d eta_m) (d eta_m / d theta): rows x q, for each class j
  DP <- lapply(seq_len(ncol(P)), function(j) Reduce(`+`, lapply(seq_along(Jm), function(m) pd$dP[, j, m] * Jm[[m]])))
  if (what == "mean_class") {
    s <- seq_len(ncol(P))
    return(list(est = W %*% (P %*% s), D = W %*% Reduce(`+`, Map(`*`, DP, s)), affected = affected,
                crossing = crossing))
  }
  list(est = W %*% P, D = do.call(rbind, lapply(DP, function(d) W %*% d)), affected = affected, crossing = crossing)
}

## ---- predictmeans() ----------------------------------------------------

## The mean type of a call and the checks of its options (one place): the
## default is "prob"; "conditional", "overall" and the options built for one
## scale stop; "mean_class" needs ordered classes (an ordinal family).
## Returns list(mean_type, parts (pm_cat_parts(), computed once: vcov() of
## a multinom fit rebuilds its Hessian from the data)).
pm_cat_check <- function(model, a) {
  pt <- pm_cat_parts(model)
  mt <- if (isTRUE(a$mean_type_default)) "prob" else a$mean_type
  if (mt %in% c("conditional", "overall")) {
    stop("mean_type = \"", mt, "\" does not apply to ", class(model)[1], " fits (a mean is a vector: one probability ",
         "per class); use mean_type = \"prob\" (the default), \"prob_diff\" or \"link\".", call. = FALSE)
  }
  if (mt == "mean_class" && !pt$ordinal) {
    stop("mean_type = \"mean_class\" needs ordered classes; a multinomial model treats them as unordered. Use ",
         "mean_type = \"prob\", or fit an ordinal model (vglm(cumulative()), clm(), polr()).", call. = FALSE)
  }
  if (mt == "link" && !is.null(pt$scale)) {
    stop("mean_type = \"link\" is not available with scale effects: the predictors are not linear in the ",
         "coefficients, so their averages over the grid are not predictors of the model. Use mean_type = \"prob\" ",
         "(the default), \"prob_diff\" or \"mean_class\".", call. = FALSE)
  }
  pm_class_opt_stop(model, "the class probabilities have their own scale (mean_type)",
                    trans = if (isTRUE(a$trans_default)) NULL else a$trans, bias_adjust = a$bias_adjust,
                    sigma = a$sigma, vcov = a$vcov, permlist = a$permlist)
  list(mean_type = mt, parts = pt)
}

## The carrier the stages of predictmeans() work on: the design of the fit
## (one coefficient per design column, for Kmatrix()'s grid and weights; the
## estimates themselves come from pm_cat_est()), the mean type and the name
## of the factor of classes (the response's name; "predictor" for "link").
pm_cat_carrier <- function(model, what, vars, pt = pm_cat_parts(model)) {
  tt <- if (inherits(model, "clm")) pm_clm_all_terms(model) else stats::terms(model)
  if (!is.null(attr(tt, "offset"))) stop("predictmeans does not handle offsets in ", class(model)[1], " fits.", call. = FALSE)
  y <- as.character(attr(tt, "variables"))[attr(tt, "response") + 1]
  var <- if (what == "link") "predictor" else if (identical(make.names(y), y)) y else "class"
  # not a name of the factors, the columns of the tables or the plots (a response 'y' would be the plot's y)
  taken <- c(vars, "y", "lo", "hi", "row", "x", "estimate", "std.error", "df", "statistic", "p.value", "conf.low",
             "conf.high", "response", "n", "estimable", "group", "by", "level1", "level2")
  if (var %in% taken) var <- paste0(var, "_class")
  structure(list(fit = model, parts = pt, what = what, var = var, call = stats::getCall(model), terms = tt,
                 contrasts = if (isS4(model)) model@contrasts else model$contrasts),
            class = "pm_catfit")
}

#' @exportS3Method stats::terms
terms.pm_catfit <- function(x, ...) x$terms
#' @exportS3Method stats::formula
formula.pm_catfit <- function(x, ...) stats::formula(x$terms)
#' @exportS3Method stats::model.frame
model.frame.pm_catfit <- function(formula, ...) stats::model.frame(formula$fit)
## the design columns only (the values are not used: pm_cat_est() replaces the means)
#' @exportS3Method NULL
pm_coef.pm_catfit <- function(model, ...) stats::setNames(rep(0, length(model$parts$xnames)), model$parts$xnames)
#' @exportS3Method NULL
pm_vcov.pm_catfit <- function(model, ...) {
  n <- model$parts$xnames
  matrix(0, length(n), length(n), dimnames = list(n, n))
}
#' @exportS3Method NULL
pm_contrasts.pm_catfit <- function(model, ...) model$contrasts
## the fit's own design columns when it has several design matrices (clm)
#' @exportS3Method NULL
pm_basis.pm_catfit <- function(model, newdata, ...) {
  if (is.null(model$parts$basis)) pm_basis.default(model, newdata, ...) else model$parts$basis(newdata)
}
## mean_type = "prob": the means are logits of the class probabilities, back-transformed as a binomial glm's
#' @exportS3Method NULL
pm_linkinv.pm_catfit <- function(model, ...) if (model$what == "prob") stats::plogis
#' @exportS3Method stats::family
family.pm_catfit <- function(object, ...) {
  if (object$what != "prob") stop("no family", call. = FALSE)
  stats::binomial()
}
#' @exportS3Method NULL
pm_df_rule.pm_catfit <- function(model, df_method, kr_info) pm_df_rule.polr(model$fit, df_method, kr_info)

## The fit's own class, wherever predictmeans() reports it.
pm_cat_fit <- function(model) if (inherits(model, "pm_catfit")) model$fit else model

## Replaces the estimates of pm_estimate() by the stacked means of the
## classes (or predictors): K the identity over (mean, class), the
## coefficients the means (logits of the probabilities for "prob"), their
## covariance D V D', z df; the factor of classes is added to the labels
## and to the data (each observation once per class: the n of a mean is its
## cell's).
pm_cat_est <- function(opt, est) {
  cf <- opt$model
  pt <- cf$parts
  what <- cf$what
  if (!is.null(est$ctr.matrix$empty_cells) || !is.null(attr(est$KK, "null_proj"))) {
    stop("Means over empty cells cannot be estimated for ", class(cf$fit)[1], " fits: the fit is rank deficient ",
         "there; leave out the interaction or the empty cells.", call. = FALSE)
  }
  K0 <- est$K
  rows <- match(rownames(K0), rownames(est$KK))
  if (what == "link") {
    ev <- pm_cat_eval(pt, K0, NULL, "link")
    lab <- pt$predictors
  } else {
    rg <- est$ctr.matrix$refgrid
    g <- as.data.frame(rg$grid)
    at <- do.call(paste, c(Map(function(nm, v) paste0(nm, " = ", if (is.numeric(v)) signif(v, 4) else v),
                               names(g), g), sep = ", "))
    ev <- pm_cat_eval(pt, rg$basis(rg$grid), as.matrix(rg$W[rows, , drop = FALSE]),
                      if (what == "mean_class") "mean_class" else "prob", labels = at)
    lab <- pt$levels
  }
  e <- as.vector(ev$est)
  D <- ev$D
  # means over a grid row where the cumulative curves cross: NA (as non-estimable means), with a caution
  aff <- if (is.null(ev$affected)) rep(FALSE, nrow(K0)) else ev$affected
  if (all(aff)) {
    stop("The fitted cumulative probabilities cross where every mean is taken: ", ev$crossing, ". The non-parallel ",
         "model is not valid there; take the means where it is (other covariate values or levels), or make the ",
         "slopes parallel.", call. = FALSE)
  }
  naff <- rep(aff, length.out = length(e))      # (stacked by class: the same means in every class)
  e[naff] <- if (what == "prob") 0.5 else 0      # (placeholders: NA in the tables, through 'estimable')
  D[naff, ] <- 0
  # separation: a class with no observations at some levels has its probability at 0 there (the MLE at
  # infinity); seen as an SE of logit(P) above 5 (that of a class observed about 0.04 times)
  n <- nrow(K0)
  pk <- pm_cat_eval(pt, K0, diag(n), "prob")
  inside <- pk$est > 0 & pk$est < 1             # (a crossing at an average design row is not separation)
  sl <- sqrt(rowSums((pk$D %*% pt$V) * pk$D)) / (pk$est * (1 - pk$est))
  sep <- which(inside & sl > 5)
  sep <- if (length(sep)) paste(rep(rownames(K0), length(pt$levels))[sep], rep(pt$levels, each = n)[sep], sep = ":")
  if (what == "prob") {                       # logit scale: d logit(P) = dP / (P (1 - P))
    if (any(e <= 0 | e >= 1)) stop("A class probability is 0 or 1: its logit is not finite; use mean_type = ",
                                   "\"prob_diff\" or \"link\".", call. = FALSE)
    D <- D / (e * (1 - e))
    e <- stats::qlogis(e)
  }
  C <- D %*% pt$V %*% t(D)
  vars <- opt$vars
  fct0 <- est$ctr.matrix$fctnames[rows, , drop = FALSE]
  if (what == "mean_class") {
    fct <- fct0
  } else {
    k <- ncol(ev$est)
    cls <- factor(rep(lab, each = n), levels = lab)
    fct <- fct0[rep(seq_len(n), k), , drop = FALSE]
    fct <- cbind(fct[vars], stats::setNames(data.frame(cls), cf$var), fct[setdiff(names(fct), vars)])
    vars <- c(vars, cf$var)
  }
  nm <- do.call(paste, c(lapply(fct[vars], as.character), sep = ":"))
  rownames(fct) <- nm
  dimnames(C) <- list(nm, nm)
  Kn <- diag(length(nm))
  dimnames(Kn) <- list(nm, nm)
  attr(Kn, "estimable") <- stats::setNames(!naff, nm)
  if (any(naff)) {
    # pairs with such a mean are not estimable either (pm_compare_sed(), as for aliased designs)
    np <- matrix(0, length(nm), sum(naff), dimnames = list(nm, NULL))
    np[cbind(which(naff), seq_len(sum(naff)))] <- 1
    attr(Kn, "null_proj") <- np
  }
  mdf <- est$mdf
  if (what != "mean_class") {
    mdf <- mdf[rep(seq_len(nrow(mdf)), k), setdiff(names(mdf), cf$var), drop = FALSE]
    mdf[[cf$var]] <- factor(rep(lab, each = nrow(est$mdf)), levels = lab)
  }
  mp <- est$mp
  mp$coef <- stats::setNames(e, nm)
  mp$vcov <- C
  mp$estimable <- rep(TRUE, length(nm))
  inf <- est$inf
  inf$df <- function(L, terms) rep(Inf, nrow(rbind(L)))
  inf$mp <- mp
  # (its message on non-estimable means does not fit here: the advice class_crossing says why)
  m <- suppressMessages(pm_estimate_means(Kn, fct, mdf, vars, mp, opt$ndecimal))
  est$ctr.matrix$K <- Kn
  est$ctr.matrix$fctnames <- fct
  # (replaced, not merged: modifyList() would merge the data frames column by column)
  new <- list(ctr.matrix = est$ctr.matrix, KK = Kn, mdf = mdf, inf = inf, mp = mp, vcovm = C, K = m$K,
              label = m$label, pm = m$pm, ses = m$ses, ndf = m$ndf, mt = m$mt, bkmt = m$mt,
              mean.table = m$mean.table, se.table = m$se.table, cat_sep = sep,
              cat_cross = if (any(aff)) list(where = ev$crossing, means = rownames(K0)[aff]))
  est[names(new)] <- new
  est
}

## The options of the later stages: the factor of classes is a factor of the
## means, and the pairs are within its levels (first, then the user's 'atvar').
pm_cat_opt <- function(opt) {
  if (opt$model$what == "mean_class") return(opt)
  v <- opt$model$var
  opt$vars <- c(opt$vars, v)
  opt$atvar <- c(v, opt$atvar)
  opt
}

## Advice of the class means (pm_cat_advice_rows(), R/advice.R) added to a
## predictmeans() result.
pm_cat_attach <- function(res, opt, est) {
  cf <- opt$model
  # the factors averaged over (covariates are held at their values, settings$covariates)
  tt <- stats::delete.response(cf$terms)
  mf <- est$mdf                                  # (the data of the means; model.frame() of a vglm can fail here)
  fac <- names(mf)[vapply(mf, function(z) is.factor(z) || is.character(z) || is.logical(z), logical(1))]
  others <- setdiff(intersect(all.vars(tt), fac), strsplit(opt$modelterm, ":", fixed = TRUE)[[1]])
  p <- if (cf$what == "prob") res$means$response else if (cf$what == "prob_diff") res$means$estimate
  out <- if (cf$what == "mean_class") any(res$means$conf.low < 1 | res$means$conf.high > length(cf$parts$levels))
  attr(res, "plot_data") <- NULL           # (the observations are classes: autoplot(data = TRUE) says so)
  # NA means from crossing curves: class_crossing says why, not the generic non-estimable line
  a <- attr(res, "advice")
  if (!is.null(est$cat_cross) && !is.null(a)) attr(res, "advice") <- a[a$id != "nonestimable_means", , drop = FALSE]
  pm_add_advice(res, pm_cat_advice_rows(cf$what, cf$parts, others, est$cat_sep, p, opt$weights, isTRUE(out),
                                        est$cat_cross))
}

## The scale line of the summary of class means ('cl': settings$response$classes).
pm_cat_scale_text <- function(cl, name) {
  switch(cl$what,
    prob = paste0("logit of the probability of each class of ", name, " (column '", cl$var, "'), back-transformed ",
                  "in 'response'; pairs within a class: log odds ratios of that class against all the others (not the ",
                  "model's coefficients, which compare classes with the reference class)"),
    prob_diff = paste0("probability of each class of ", name, " (column '", cl$var, "'); pairs within a class: ",
                       "differences of probabilities (for two means they sum to 0 over the classes)"),
    link = paste0("the linear predictors of the model (column '", cl$var, "'); pairs within a predictor"),
    mean_class = paste0("mean class of ", name, " (the classes scored 1, 2, 3, ...)"))
}

## The stop of the functions that work on one scale (everything but
## predictmeans()): reached through pm_coef() and the ordinal checks.
pm_cat_stop <- function(model) {
  stop(class(model)[1], " fits of a class response without a single latent scale (multinomial, slopes that are ",
       "not parallel, or clm nominal or scale effects) have a vector of means: predictmeans() gives their class probabilities (mean_type = \"prob\"); ",
       "this function does not handle them.", call. = FALSE)
}

#' @exportS3Method NULL
pm_coef.multinom <- function(model, ...) pm_cat_stop(model)
#' @exportS3Method NULL
pm_vcov.multinom <- function(model, ...) pm_cat_stop(model)

## ---- clm with nominal or scale effects (B4 G1) -----------------------

## P(Y <= j) = F((theta_j(z) - x'b) / exp(s'zeta)): theta_j(z) = sum_c z_c
## (tJac alpha_c)_j, alpha_c the threshold parameters of nominal column c
## (alpha.mat; tJac maps them to the J - 1 cut points). In theta (coef()
## order: alpha by nominal column, then b, then zeta) the numerator is
## linear: H[, alpha_cl] = tJac[, l] on column c, H[, b] = -1; zeta enters
## through 'scale'. The design columns: the location columns of b, then
## "nominal:" and "scale:" columns (a variable can be in several formulas:
## clm names its coefficients alike, so they are matched by position).
pm_clm_cat_parts <- function(model) {
  pm_clm_link_check(model)
  al <- unlist(model$aliased)
  if (any(al)) {
    stop("The clm fit has aliased coefficients (", paste(names(al)[al], collapse = ", "), "): the means would not be ",
         "unique. Refit without the aliased terms (a variable in both 'formula' and 'nominal' belongs in one).",
         call. = FALSE)
  }
  am <- if (is.null(model$alpha.mat)) {                  # nominal columns x threshold parameters (none: the intercept)
    matrix(model$alpha, 1, dimnames = list("(Intercept)", names(model$alpha)))
  } else as.matrix(model$alpha.mat)
  tJ <- as.matrix(model$tJac)                           # cut points x threshold parameters
  b <- model$beta
  z <- if (length(model$zeta)) model$zeta else numeric(0)
  theta <- c(as.vector(t(am)), b, z)
  V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
  cf <- stats::coef(model)
  if (length(cf) != length(theta) || max(abs(unname(cf) - theta)) > 1e-8 * max(1, abs(theta))) {
    stop("internal: the clm coefficients do not match alpha.mat, beta and zeta.", call. = FALSE)
  }
  if (is.null(V) || anyNA(V) || any(!is.finite(V)) || nrow(V) != length(theta)) {
    stop("The clm fit has no usable covariance matrix of its coefficients (refit with Hess = TRUE; a class that ",
         "never occurs in some cells?).", call. = FALSE)
  }
  names(theta) <- rownames(V) <- colnames(V) <- make.unique(names(cf))
  nn <- paste0("nominal:", rownames(am))
  sn <- if (length(z)) paste0("scale:", names(z)) else character(0)
  xnames <- c(names(b), nn, sn)
  M <- nrow(tJ)
  L <- ncol(am)
  H <- cbind(do.call(cbind, rep(list(tJ), nrow(am))), matrix(-1, M, length(b)), matrix(0, M, length(z)))
  col <- c(rep(length(b) + seq_len(nrow(am)), each = L), seq_along(b), length(b) + nrow(am) + seq_along(z))
  mm <- function(tt, ctr, xlev, nd, keep) {
    tt <- stats::delete.response(tt)
    m <- stats::model.frame(tt, nd, na.action = stats::na.pass, xlev = xlev)
    stats::model.matrix(tt, m, contrasts.arg = ctr)[, keep, drop = FALSE]
  }
  basis <- function(nd) {
    B <- cbind(mm(model$terms, model$contrasts, model$xlevels, nd, names(b)),
               if (is.null(model$nom.terms)) matrix(1, nrow(nd), 1)
               else mm(model$nom.terms, model$nom.contrasts, model$nom.xlevels, nd, rownames(am)),
               if (length(z)) mm(model$S.terms, model$S.contrasts, model$S.xlevels, nd, names(z)))
    colnames(B) <- xnames
    B
  }
  Fd <- pm_ord_link_fd(model$link)
  J <- M + 1
  pd <- function(eta) {
    Fe <- Fd[[1]](eta)
    fe <- Fd[[2]](eta)
    P <- cbind(Fe, 1) - cbind(0, Fe)
    dP <- array(0, c(nrow(P), J, M))
    for (m in seq_len(M)) {
      dP[, m, m] <- fe[, m]
      dP[, m + 1, m] <- -fe[, m]
    }
    list(P = P, dP = dP)
  }
  lev <- model$y.levels
  list(theta = theta, V = V, H = H, col = col, xnames = xnames, pd = pd, levels = lev,
       predictors = paste0(model$link, "(P(Y <= ", lev[-J], "))"),
       reference = NULL, family = "cumulative", ordinal = TRUE, basis = basis,
       scale = if (length(z)) list(cols = length(b) + nrow(am) + seq_along(z), idx = length(theta) - length(z) + seq_along(z)))
}

## The terms of all the variables of a clm fit (location, nominal, scale),
## for the grid of the means.
pm_clm_all_terms <- function(model) {
  parts <- list(model$terms, model$nom.terms, model$S.terms)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (any(vapply(parts, function(t) !is.null(attr(t, "offset")), logical(1)))) {
    stop("predictmeans does not handle offsets in clm fits.", call. = FALSE)
  }
  labs <- unique(unlist(lapply(parts, attr, "term.labels")))
  y <- as.character(attr(model$terms, "variables"))[attr(model$terms, "response") + 1]
  f <- stats::reformulate(if (length(labs)) labs else "1", response = y)
  environment(f) <- environment(stats::formula(model$terms))
  stats::terms(f)
}
