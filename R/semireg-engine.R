## The fitting engine of semireg() and semireg_tmb() (dev/review-semireg.md).
##
## A semi-parametric model is a mixed model in which each smooth is a block
## of random effects with a given design Z (from smZ()) and iid coefficients.
## Both front ends:
##   1. pm_sr_setup(): the data (complete cases, numeric covariates centred
##      or scaled), the smooth designs, and the formula with one placeholder
##      term (1 | name) per smooth block;
##   2. build the model structure WITHOUT fitting (lme4::lFormula() or
##      glFormula(); glmmTMB(doFit = FALSE)), replace the placeholder rows of
##      Z by the smooth designs (pm_sr_replace_Zt()), and fit once;
##   3. pm_sr_info(): the covariance of the fixed and smooth coefficients and
##      what semipred() needs to predict.
## The data are given to the builders in an environment of their own, never
## written into the caller's frame.

## The data, the smooth designs and the formula with placeholders.
##   mc         the front end's match.call() (formula already set)
##   smoothZ    the unevaluated smoothZ argument
##   env        where the front end was called
## Returns list(data, fomul_vars, sm_vars, smoothZm (as given, by name),
## smoothZN (blocks: a smZ(group = TRUE) list is split into its levels),
## Zt (transposed blocks), formula (with placeholders, or NULL without
## smooths)).
pm_sr_setup <- function(mc, data, smoothZ, env, ncenter = FALSE, nscale = FALSE, resp_scale = FALSE) {
  formula <- mc$formula
  # smooths written in the formula (R/semireg-smooth.R)
  sm_f <- pm_sm_parse(formula, data)
  fomul_vars <- all.vars(stats::terms(if (is.null(sm_f)) formula else sm_f$formula))
  sm_vars <- intersect(all.vars(smoothZ), names(data))
  if (any(unlist(lapply(data[, setdiff(sm_vars, fomul_vars), drop = FALSE], is.factor)))) {
    stop("Any factor in 'smoothZ' list must be in the formula!")
  }
  # the response may have several variables (cbind(successes, failures)):
  # none of them is centred or scaled unless resp_scale = TRUE
  response_vars <- all.vars(formula[[2]])
  fomul_vars <- fomul_vars[-1]
  data <- droplevels(as.data.frame(stats::na.omit(data[, intersect(all.vars(mc), names(data)), drop = FALSE])))
  numeric_n <- names(data)[vapply(data, is.numeric, logical(1))]
  if (!resp_scale) numeric_n <- setdiff(numeric_n, response_vars)
  if (!is.null(sm_f)) numeric_n <- setdiff(numeric_n, sm_f$vars)   # variables of formula smooths: as they are
  if (any(ncenter, nscale)) {
    for (i in numeric_n) {
      if (!is.matrix(data[, i])) {
        scaled_i <- scale(data[i], center = ncenter, scale = nscale)
        data[i] <- as.numeric(scaled_i)
        if (ncenter) attr(data, paste(i, "mean", sep = "_")) <- attr(scaled_i, "scaled:center")
        if (nscale) attr(data, paste(i, "sd", sep = "_")) <- attr(scaled_i, "scaled:scale")
      }
    }
    attr(data, "numeric_var") <- numeric_n
  }
  smoothZm <- eval(smoothZ, data, env)
  smooths <- NULL
  if (!is.null(sm_f) && length(smoothZm)) {
    stop("Give the smooths either in the formula (s(), t2()) or in 'smoothZ', not both.", call. = FALSE)
  }
  if (!is.null(sm_f)) {
    # the smooths of the formula: random blocks, and their fixed terms added
    # to the formula (those already there are not doubled)
    b <- pm_sm_build(sm_f, data)
    smoothZm <- c(smoothZm, b$blocks)
    smooths <- b$smooths
    sm_vars <- union(sm_vars, b$vars)
    lab <- attr(stats::terms(sm_f$formula), "term.labels")
    add <- setdiff(b$fixed, lab)
    formula <- if (length(add)) {
      stats::as.formula(paste(paste(deparse(sm_f$formula, width.cutoff = 500L), collapse = " "), "+",
                              paste(add, collapse = " + ")), env = environment(formula))
    } else sm_f$formula
    fomul_vars <- all.vars(stats::terms(formula))[-1]
  }
  smoothZN <- list()
  for (i in names(smoothZm)) {
    smoothZN <- c(smoothZN, if (is.list(smoothZm[[i]])) unlist(smoothZm[i], recursive = FALSE) else smoothZm[i])
  }
  Zt <- lapply(smoothZN, t)
  if (length(Zt)) {
    stopifnot(length(names(Zt)) == length(Zt), all(vapply(Zt, methods::is, logical(1), class2 = "sparseMatrix")))
  }
  # one placeholder random intercept per smooth block, its factor in the data
  f_new <- NULL
  if (length(Zt)) {
    rhs <- paste(paste(format(formula), collapse = ""), paste0("+ (1|", names(Zt), ")", collapse = " "))
    for (i in names(Zt)) data[[i]] <- factor(rep(seq_len(nrow(Zt[[i]])), length = ncol(Zt[[i]])))
    f_new <- stats::as.formula(rhs, env = environment(formula))
  }
  list(data = data, fomul_vars = fomul_vars, sm_vars = sm_vars, smoothZm = smoothZm, smoothZN = smoothZN, Zt = Zt,
       formula = f_new, smooths = smooths)
}

## Evaluate a builder call with the semireg data, in an environment of its
## own whose parent is the caller's (for weights, offsets, subsets that refer
## to other objects).
pm_sr_eval <- function(cl, data, env) {
  e <- new.env(parent = env)
  e$.pm_sr_data <- data
  cl$data <- quote(.pm_sr_data)
  eval(cl, e)
}

## Replace the placeholder rows of Zt (one random-intercept term per smooth
## block, found through flist/Gp) by the smooth designs.
pm_sr_replace_Zt <- function(Zt, smoothZt, flist, Gp) {
  asgn <- attr(flist, "assign")
  stopifnot(all(names(smoothZt) %in% names(flist)))
  for (i in names(smoothZt)) {
    tn <- which(match(i, names(flist)) == asgn)
    if (length(tn) > 1) stop("a smoothZt factor must be associated with only one r.e. term")
    rowsi <- (Gp[tn] + 1L):Gp[tn + 1L]
    stopifnot(all(dim(Zt[rowsi, ]) == dim(smoothZt[[i]])))
    Zt[rowsi, ] <- smoothZt[[i]]
  }
  Zt
}

## Keep only the arguments of a call that the function 'fn' accepts by name.
pm_sr_call_for <- function(cl, fn, drop = character(0)) {
  keep <- c("", intersect(names(cl)[-1], setdiff(names(formals(fn)), c("...", drop))))
  cl[names(cl) %in% keep | seq_along(cl) == 1L]
}

## lme4: the structure (lFormula()/glFormula()), the smooth designs put into
## Zt, one fit. 'lmerc' is the [g]lmer call of the front end (formula with
## placeholders). Returns the merMod and, for a gaussian fit, the deviance
## function that goes with it (lmerTest needs it for Satterthwaite df).
pm_sr_fit_lme4 <- function(lmerc, st, gaus, family, ngbinomial, control, start, verbose, drop.unused.levels, mc, env,
                           ...) {
  fc <- lmerc
  fc$formula <- st$formula
  if (gaus) {
    fc[[1]] <- quote(lme4::lFormula)
    fc <- pm_sr_call_for(fc, lme4::lFormula)
  } else {
    fc[[1]] <- quote(lme4::glFormula)
    fc$family <- if (ngbinomial) stats::poisson() else family
    fc <- pm_sr_call_for(fc, lme4::glFormula, drop = "start")
  }
  lf <- suppressWarnings(pm_sr_eval(fc, st$data, env))
  reTrms <- lf$reTrms
  reTrms$Zt <- pm_sr_replace_Zt(reTrms$Zt, st$Zt, reTrms$flist, reTrms$Gp)
  dfl <- list(fr = lf$fr, X = lf$X, reTrms = reTrms, start = reTrms$theta, drop.unused.levels = drop.unused.levels)
  glmm_fit <- function(fam) {
    dfl$family <- fam
    devfun <- do.call(lme4::mkGlmerDevfun, dfl)
    opt <- lme4::optimizeGlmer(devfun, start = start, verbose = verbose, optimizer = control$optimizer,
                               control = control$optCtrl, ...)
    devfun <- lme4::updateGlmerDevfun(devfun, dfl$reTrms)
    opt <- lme4::optimizeGlmer(devfun, stage = 2, start = start, verbose = verbose, optimizer = control$optimizer,
                               control = control$optCtrl, ...)
    lme4::mkMerMod(environment(devfun), opt, reTrms, lf$fr, mc)
  }
  if (ngbinomial) {
    # negative binomial: theta by ML from a Poisson fit, then the fit with it
    semer1 <- glmm_fit(stats::poisson())
    Y <- stats::model.response(stats::model.frame(semer1))
    mu <- stats::na.omit(stats::fitted(semer1))
    theta_v <- MASS::theta.ml(Y, mu, weights = semer1@resp$weights, limit = 20)
    return(list(semer = glmm_fit(MASS::negative.binomial(theta = theta_v)), devfun = NULL))
  }
  if (!gaus) return(list(semer = glmm_fit(family), devfun = NULL))
  dfl$REML <- lf$REML
  devfun <- do.call(lme4::mkLmerDevfun, dfl)
  opt <- lme4::optimizeLmer(devfun, start = start, verbose = verbose, optimizer = control$optimizer,
                            control = control$optCtrl, ...)
  list(semer = lme4::mkMerMod(environment(devfun), opt, reTrms, lf$fr, mc), devfun = devfun)
}

## glmmTMB: the structure (doFit = FALSE), the smooth designs put into Z,
## one fit.
pm_sr_fit_tmb <- function(tmbc, st, env) {
  tmbc$formula <- st$formula
  tmbc$doFit <- FALSE
  tmbf <- pm_sr_eval(tmbc, st$data, env)
  rt <- tmbf$condList$reTrms
  tmbf$data.tmb$Z <- t(pm_sr_replace_Zt(rt$Zt, st$Zt, rt$flist, rt$Gp))
  list(semer = glmmTMB::fitTMB(tmbf), tmbf = tmbf)
}

## The covariance of the fixed and smooth coefficients, relative to sigma^2,
## and what semipred() needs (the coefficients of each smooth, its knots and
## type). With b = Lambda u, u ~ N(0, sigma^2 I) (as lme4), C~ = [X, Z Lambda]
## and T = diag(I, Lambda):
##   Cov(beta, b)/sigma^2 = T (C~' W C~ + diag(0, I))^-1 T'
## which equals (C' W C + D)^-1 with D = (Lambda Lambda')^-1 when Lambda is
## invertible, and needs no inverse of Lambda: a smooth whose variance is
## estimated as zero (a straight line) gets zero covariance. (Before 2.0 the
## second form was used, with nearPD() for a singular Lambda, which failed.)
##   X, Z, Lambda   of the fit; w: weights (length n)
##   uHat           random effects in the column order of Z, named
##   u_ind          number of columns of Z of each term, named by term
pm_sr_info <- function(X, Z, Lambda, w, uHat, u_ind, Zt_name, smoothZm) {
  n_beta <- ncol(X)
  C <- cbind(X, Z)
  CTC <- Matrix::crossprod(C * w, C)
  Ct <- cbind(X, Z %*% Lambda)
  q <- ncol(Z)
  Tm <- Matrix::.bdiag(list(Matrix::Diagonal(n_beta), Lambda))
  M <- Matrix::crossprod(Ct * w, Ct) + Matrix::.bdiag(list(Matrix::Matrix(0, n_beta, n_beta), Matrix::Diagonal(q)))
  fullCovMat <- Tm %*% Matrix::solve(M, Matrix::t(Tm))
  fullCovMat <- Matrix::forceSymmetric(methods::as(fullCovMat, "generalMatrix"))
  df <- sum(Matrix::diag(fullCovMat %*% CTC))              # effective degrees of freedom

  vcov_ind <- cumsum(c(n_beta, u_ind))
  names(vcov_ind) <- c("beta", names(u_ind))
  Zt_name <- intersect(names(u_ind), Zt_name)
  uHat_ind <- vcov_ind - n_beta
  cov_lst <- u_lst <- vector("list", length(Zt_name))
  names(cov_lst) <- names(u_lst) <- Zt_name
  block <- function(nm, ind) (ind[match(nm, names(ind)) - 1] + 1):ind[match(nm, names(ind))]
  for (i in Zt_name) {
    cov_lst[[i]] <- fullCovMat[c(seq_len(n_beta), block(i, vcov_ind)), c(seq_len(n_beta), block(i, vcov_ind))]
    u_lst[[i]] <- uHat[block(i, uHat_ind)]
  }
  # a smZ(group = TRUE) term: its blocks together
  vcov_indN <- vcov_ind
  names(vcov_indN) <- gsub("\\..*", "", names(vcov_ind))
  vcov_indN <- rev(vcov_indN)[unique(names(vcov_indN))]
  for (i in setdiff(names(vcov_indN), names(vcov_ind))) {
    idx <- c(seq_len(n_beta), block(i, vcov_indN))
    cov_lst[[i]] <- fullCovMat[idx, idx]
  }
  sm_idx <- c(seq_len(n_beta), unlist(lapply(Zt_name, block, ind = vcov_ind)))
  cov_lst$sm_Cov <- fullCovMat[sm_idx, sm_idx]

  Z_namem <- intersect(names(vcov_indN), names(smoothZm))
  knots_lst <- range_lst <- type_lst <- stats::setNames(vector("list", length(Z_namem)), Z_namem)
  for (i in Z_namem) {
    knots_lst[[i]] <- attr(smoothZm[[i]], "knots")
    range_lst[[i]] <- attr(smoothZm[[i]], "range.x")
    type_lst[[i]] <- attr(smoothZm[[i]], "type")
  }
  list(knots_lst = knots_lst, range_lst = range_lst, type_lst = type_lst, cov_lst = cov_lst, u_lst = u_lst,
       CovMat = fullCovMat, Cov_ind = vcov_ind, Cov_indN = vcov_indN, df = df)
}

## The "semireg" result of both front ends.
pm_sr_result <- function(semer, st, mc, info, extra) {
  if (is.null(info)) info <- list(knots_lst = NULL, range_lst = NULL, type_lst = NULL, cov_lst = NULL, u_lst = NULL,
                                  CovMat = NULL, Cov_ind = NULL, Cov_indN = NULL, df = NULL)
  ans <- c(list(semer = semer, data = st$data, fomul_vars = st$fomul_vars, sm_vars = st$sm_vars,
                smoothZ_call = mc$smoothZ), info, extra, list(smooths = st$smooths))
  class(ans) <- c("semireg", "list")
  ans
}
