## Categorical responses without a single latent scale (B3b, dev/review-vglm.md
## E2): nnet::multinom(), VGAM::vglm() with multinomial(), and vglm ordinal
## families (cumulative, acat, cratio, sratio) whose slopes are not parallel.
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
## coefficient). probs(eta): the class probabilities; deta(eta): their
## Jacobian in eta (J x M). levels: the classes; predictors: the names of
## the linear predictors; reference: the reference class of a multinomial
## fit; ordinal: TRUE for the ordinal families.
pm_cat_parts <- function(model) {
  if (inherits(model, "multinom")) {
    if (isTRUE(model$censored) || !isTRUE(model$softmax) && !isTRUE(model$entropy)) {
      stop("predictmeans handles multinom fits of a class response (softmax, or a two-class response); not censored ",
           "or other nnet fits.", call. = FALSE)
    }
    if (isTRUE(model$decay > 0)) {
      stop("The multinom fit is penalised (decay = ", model$decay, "): its Hessian is not the information of the ",
           "estimates, so SEs and tests would be wrong; refit with decay = 0.", call. = FALSE)
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
    probs <- function(eta) { e <- exp(c(0, eta) - max(0, eta)); e / sum(e) }
    deta <- function(eta) { P <- probs(eta); (diag(P, length(P)) - tcrossprod(P))[, -1, drop = FALSE] }
    ref <- lev[1]
    predictors <- paste0("log(P[", lev[-1], "]/P[", ref, "])")
    family <- "multinomial"
  } else {
    fam <- model@family@vfamily[1]
    if (!is.null(model@control$xij)) {
      stop("vglm fits with 'xij' (covariates that differ between predictors) are not supported.", call. = FALSE)
    }
    if (any(grepl("^sm\\.", all.names(stats::formula(model))))) {
      stop("VGAM's smart terms (sm.bs, sm.ns, sm.poly, sm.scale) cannot be rebuilt on the grid of the means; refit ",
           "with bs(), ns(), poly() or scale().", call. = FALSE)
    }
    cm <- VGAM::constraints(model)
    theta <- VGAM::coef(model)
    V <- tryCatch(VGAM::vcov(model), error = function(e) NULL)
    H <- do.call(cbind, unname(cm))
    col <- rep(seq_along(cm), vapply(cm, ncol, integer(1)))
    xnames <- names(cm)
    if (ncol(H) != length(theta)) stop("internal: the vglm constraints do not match its coefficients.", call. = FALSE)
    # (two rows: some families' linkinv drop a one-row matrix)
    probs <- function(eta) as.vector(model@family@linkinv(rbind(eta, eta), extra = model@extra)[1, ])
    deta <- function(eta) pm_fd_jac(probs, eta)
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
  list(theta = theta, V = V, H = H, col = col, xnames = xnames, probs = probs, deta = deta, levels = lev,
       predictors = predictors, reference = ref, family = family, ordinal = family != "multinomial")
}

## Means and their gradient in theta. X: rows of the design (columns
## xnames). what = "link": the predictors at the rows of X (linear: X are
## the K rows); "prob" / "mean_class": at the grid rows X, averaged with the
## weights W (rows: the means). Returns list(est, D): est a matrix (means x
## classes, predictors, or one column), D the gradient of as.vector(est)
## (stacked by column: all means of the first class, then the second, ...).
pm_cat_eval <- function(pt, X, W = NULL, what = c("prob", "link", "mean_class")) {
  what <- match.arg(what)
  X <- as.matrix(X)[, pt$xnames, drop = FALSE]
  Xc <- X[, pt$col, drop = FALSE]                          # x[col] for each coefficient
  M <- nrow(pt$H)
  if (what == "link") {
    est <- Xc %*% (pt$theta * t(pt$H))
    D <- do.call(rbind, lapply(seq_len(M), function(j) sweep(Xc, 2, pt$H[j, ], "*")))
    return(list(est = est, D = D))
  }
  eta <- Xc %*% (pt$theta * t(pt$H))                        # grid rows x M
  J <- length(pt$levels)
  P <- matrix(0, nrow(X), J)
  DP <- array(0, c(nrow(X), J, length(pt$theta)))
  for (g in seq_len(nrow(X))) {
    P[g, ] <- pt$probs(eta[g, ])
    # d P / d theta = d P / d eta  x  d eta / d theta (= H with each column times x[col])
    DP[g, , ] <- pt$deta(eta[g, ]) %*% sweep(pt$H, 2, Xc[g, ], "*")
  }
  # non-parallel cumulative fits: the cumulative curves can cross, giving negative "probabilities"
  bad <- which(P < -1e-10, arr.ind = TRUE)
  if (nrow(bad)) {
    stop("The fitted cumulative probabilities cross where the means are taken (class '", pt$levels[bad[1, 2]],
         "' has probability ", format(signif(P[bad[1, ]], 3)), " at row ", bad[1, 1], " of the grid): the ",
         "non-parallel model is not valid there. Take the means within the range where it is valid (factor levels ",
         "or covariate values), or make the slopes parallel.", call. = FALSE)
  }
  if (what == "mean_class") {
    s <- seq_len(J)
    m <- P %*% s
    Dm <- apply(DP, 3, function(d) d %*% s)
    return(list(est = W %*% m, D = W %*% matrix(Dm, nrow(X))))
  }
  list(est = W %*% P, D = do.call(rbind, lapply(seq_len(J), function(j) W %*% matrix(DP[, j, ], nrow(X)))))
}

## ---- predictmeans() ----------------------------------------------------

## The mean type of a call and the checks of its options (one place): the
## default is "prob"; "conditional", "overall" and the options built for one
## scale stop; "mean_class" needs ordered classes (an ordinal family).
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
  bad <- c(trans = !isTRUE(a$trans_default) && !is.null(a$trans) && !identical(a$trans, I),
           bias_adjust = !isFALSE(a$bias_adjust), sigma = !is.null(a$sigma), vcov = !is.null(a$vcov),
           permlist = !is.null(a$permlist))
  if (any(bad)) {
    stop("'", paste(names(bad)[bad], collapse = "', '"), "' cannot be used with ", class(model)[1], " fits: the ",
         "class probabilities have their own scale (mean_type).", call. = FALSE)
  }
  mt
}

## The carrier the stages of predictmeans() work on: the design of the fit
## (one coefficient per design column, for Kmatrix()'s grid and weights; the
## estimates themselves come from pm_cat_est()), the mean type and the name
## of the factor of classes (the response's name; "predictor" for "link").
pm_cat_carrier <- function(model, what, vars) {
  pt <- pm_cat_parts(model)
  tt <- stats::terms(model)
  if (!is.null(attr(tt, "offset"))) stop("predictmeans does not handle offsets in ", class(model)[1], " fits.", call. = FALSE)
  y <- as.character(attr(tt, "variables"))[attr(tt, "response") + 1]
  var <- if (what == "link") "predictor" else if (identical(make.names(y), y)) y else "class"
  if (var %in% vars) var <- paste0(".", var)
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
    ev <- pm_cat_eval(pt, rg$basis(rg$grid), as.matrix(rg$W[rows, , drop = FALSE]),
                      if (what == "mean_class") "mean_class" else "prob")
    lab <- pt$levels
  }
  e <- as.vector(ev$est)
  D <- ev$D
  # separation: a class with no observations at some levels has its probability at 0 there (the MLE at
  # infinity); seen as an SE of logit(P) above 5 (that of a class observed about 0.04 times)
  n <- nrow(K0)
  pk <- pm_cat_eval(pt, K0, diag(n), "prob")
  sl <- sqrt(rowSums((pk$D %*% pt$V) * pk$D)) / (pk$est * (1 - pk$est))
  sep <- which(!is.finite(sl) | sl > 5)
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
  attr(Kn, "estimable") <- stats::setNames(rep(TRUE, length(nm)), nm)
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
  m <- pm_estimate_means(Kn, fct, mdf, vars, mp, opt$ndecimal)
  est$ctr.matrix$K <- Kn
  est$ctr.matrix$fctnames <- fct
  # (replaced, not merged: modifyList() would merge the data frames column by column)
  new <- list(ctr.matrix = est$ctr.matrix, KK = Kn, mdf = mdf, inf = inf, mp = mp, vcovm = C, K = m$K,
              label = m$label, pm = m$pm, ses = m$ses, ndf = m$ndf, mt = m$mt, bkmt = m$mt,
              mean.table = m$mean.table, se.table = m$se.table, cat_sep = sep)
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
  tt <- stats::delete.response(stats::terms(cf$fit))
  others <- setdiff(all.vars(tt), strsplit(opt$modelterm, ":", fixed = TRUE)[[1]])
  p <- if (cf$what == "prob") res$means$response else if (cf$what == "prob_diff") res$means$estimate
  attr(res, "plot_data") <- NULL           # (the observations are classes: autoplot(data = TRUE) says so)
  pm_add_advice(res, pm_cat_advice_rows(cf$what, cf$parts, others, est$cat_sep, p))
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
  stop(class(model)[1], " fits of a class response without a single latent scale (multinomial, or slopes that are ",
       "not parallel) have a vector of means: predictmeans() gives their class probabilities (mean_type = \"prob\"); ",
       "this function does not handle them.", call. = FALSE)
}

#' @exportS3Method NULL
pm_coef.multinom <- function(model, ...) pm_cat_stop(model)
#' @exportS3Method NULL
pm_vcov.multinom <- function(model, ...) pm_cat_stop(model)
