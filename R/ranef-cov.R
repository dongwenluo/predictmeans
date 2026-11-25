## The random effects of a mixed model and the joint covariance of
## (beta_hat, u_hat - u): the inverse of the mixed model equations (Henderson's
## C; SAS "Std Err Pred"), on the absolute scale (sigma^2 included), for
## se_ranef(). One place for every class:
##   lmer, glmer   pm_sr_info() on X, Z, Lambda and the (working) weights
##   glmmTMB       the (beta, b) block of TMB's joint precision, inverted
##                 (pm_tmb_cond_cov(); given the variance parameters)
##   lme           V = Z G Z' + R rebuilt by pm_nlme_reml() (any reStruct,
##                 corStruct, varFunc; checked against vcov()), then
##                 C11 = (X'V^-1X)^-1, C12 = -C11 X'V^-1 Z G,
##                 C22 = G - G Z'V^-1 Z G + G Z'V^-1 X C11 X'V^-1 Z G
##   semireg       its fit$semer (one of the above)
## cond: the conditional variances Var(u | y), beta and the variance
## parameters known (lme4 condVar; for lme G - G Z'V^-1 Z G); for glmmTMB
## the diagonal of the inverse of the whole joint precision (what
## ranef(condVar = TRUE) gives: beta and variance parameters NOT known).
## Returns list(beta, u, V (p + q square), cond (q), tab (group, level, term,
## variance: the term's variance component), rows (per group: the factor of
## the data rows), data, fixed (fixed-effect factors), family).
## (One method per class: merMod, glmmTMB, lme; semireg passes on its fit.)
pm_ranef_cov <- function(model) UseMethod("pm_ranef_cov")
#' @exportS3Method NULL
pm_ranef_cov.semireg <- function(model) pm_ranef_cov(pm_fit(model))
#' @exportS3Method NULL
pm_ranef_cov.aovlist <- function(model) {
  stop("se_ranef() needs a model with random effects; an aov with Error() has strata, not predicted random ",
       "effects: fit the strata as random effects (lmer(), lme()).", call. = FALSE)
}
#' @exportS3Method NULL
pm_ranef_cov.default <- function(model) {
  stop("se_ranef() handles models with random effects fitted by lmer, glmer, glmmTMB, lme or semireg; a ",
       class(model)[1], " model is not covered.", call. = FALSE)
}

## Rows of the table of random effects for lme4 and glmmTMB fits: Z columns run
## by term, within a term by level, then component.
pm_ranef_tab_mer <- function(model, flist, re) {
  asg <- attr(flist, "assign")
  do.call(rbind, lapply(seq_along(re), function(j) {
    tm <- re[[j]]
    lev <- levels(flist[[asg[j]]])
    nc <- length(tm$comps)
    data.frame(group = names(flist)[asg[j]], level = rep(lev, each = nc), term = rep(tm$comps, length(lev)),
               variance = rep(diag(tm$vcov), length(lev)), stringsAsFactors = FALSE)
  }))
}

#' @exportS3Method NULL
pm_ranef_cov.merMod <- function(model) {
  X <- lme4::getME(model, "X")
  Z <- lme4::getME(model, "Z")
  Lambda <- lme4::getME(model, "Lambda")
  b <- lme4::getME(model, "b")
  s2 <- stats::sigma(model)^2
  gaus <- lme4::isLMM(model)
  w <- if (gaus) model@resp$weights else model@resp$sqrtWrkWt()^2
  # lme4's GLMM working weights include 1/sigma^2 (pm_sr_info(s2 = ))
  info <- pm_sr_info(X, Z, Lambda, w, b, c(all = ncol(Z)), character(0), s2 = if (gaus) 1 else s2)
  V <- as.matrix(info$CovMat) * s2
  # conditional (lme4 condVar): sigma^2 Lambda (Lambda'Z'WZ Lambda + I)^-1 Lambda'
  ZL <- Z %*% Lambda
  Mq <- Matrix::forceSymmetric(Matrix::crossprod(ZL * w, ZL) + Matrix::Diagonal(ncol(Z)))
  cond <- Matrix::diag(Lambda %*% Matrix::solve(Mq, Matrix::t(Lambda))) * (if (gaus) s2 else 1)
  flist <- lme4::getME(model, "flist")
  re <- pm_re_terms(model)
  tab <- pm_ranef_tab_mer(model, flist, re)
  mf <- stats::model.frame(model)
  list(beta = lme4::fixef(model)[colnames(X)], u = as.vector(b), V = V, cond = as.vector(cond), tab = tab,
       rows = as.list(flist)[unique(names(flist)[attr(flist, "assign")])], data = mf,
       fixed = pm_ranef_fixed_factors(stats::formula(model, fixed.only = TRUE), mf),
       family = stats::family(model), cond_label = "lme4 condVar: Var(u | y), fixed effects known")
}

#' @exportS3Method NULL
pm_ranef_cov.glmmTMB <- function(model) {
  pm_need("TMB", "se_ranef() of glmmTMB fits")
  X <- lme4::getME(model, "X")
  q <- ncol(lme4::getME(model, "Z"))
  P <- TMB::sdreport(model$obj, getJointPrecision = TRUE)$jointPrecision
  nm <- rownames(P)
  ib <- which(nm == "beta")
  iu <- which(nm == "b")[seq_len(q)]
  idx <- c(ib, iu)
  V <- as.matrix(solve(P[idx, idx]))
  cond <- diag(as.matrix(solve(P)))[iu]
  flist <- model$modelInfo$reTrms$cond$flist
  re <- pm_re_terms(model)
  tab <- pm_ranef_tab_mer(model, flist, re)
  mf <- stats::model.frame(model)
  list(beta = glmmTMB::fixef(model)$cond[colnames(X)], u = as.vector(lme4::getME(model, "b")[seq_len(q)]), V = V,
       cond = cond, tab = tab, rows = as.list(flist)[unique(names(flist)[attr(flist, "assign")])], data = mf,
       fixed = pm_ranef_fixed_factors(stats::formula(model, fixed.only = TRUE), mf),
       family = stats::family(model),
       cond_label = "glmmTMB condVar (TMB sdreport): includes the uncertainty of the fixed effects and variance parameters")
}

#' @exportS3Method NULL
pm_ranef_cov.lme <- function(model) {
  R <- pm_nlme_reml(model)                     # stops when V cannot be rebuilt
  Srel <- R$Srel
  Zs <- attr(Srel, "Zs")
  th <- R$th0[-R$k]
  D <- attr(Srel, "pieces")(th, "D")$D
  s2 <- model$sigma^2
  dat <- pm_get_data(model)
  grps <- nlme::getGroups(dat, nlme::getGroupsFormula(model))
  if (!is.data.frame(grps)) {
    grps <- data.frame(grps)
    names(grps) <- names(Zs)[1]
  }
  re_f <- stats::formula(model$modelStruct$reStruct)
  rows <- list()
  tab <- list()
  for (lv in names(Zs)) {
    gq <- interaction(grps[seq_len(match(lv, names(grps)))], drop = TRUE, sep = "/")
    rows[[lv]] <- gq
    Dl <- as.matrix(D[[lv]])
    cm <- colnames(stats::model.matrix(re_f[[lv]], dat))
    tab[[lv]] <- data.frame(group = lv, level = rep(levels(gq), each = length(cm)), term = rep(cm, nlevels(gq)),
                            variance = rep(diag(Dl) * s2, nlevels(gq)), stringsAsFactors = FALSE)
  }
  Z <- do.call(cbind, lapply(Zs, `[[`, "Z"))
  G <- Matrix::bdiag(lapply(names(Zs), function(lv) Matrix::kronecker(Zs[[lv]]$I, Matrix::Matrix(D[[lv]]))))
  X <- R$X
  Ch <- Matrix::Cholesky(Srel(th), LDL = FALSE)
  ViX <- as.matrix(Matrix::solve(Ch, X, system = "A"))
  ZG <- Z %*% G
  ViZG <- Matrix::solve(Ch, ZG, system = "A")
  C11 <- solve(crossprod(X, ViX))
  XViZG <- as.matrix(Matrix::crossprod(X, ViZG))          # p x q
  C12 <- -C11 %*% XViZG
  Cc <- as.matrix(G - Matrix::crossprod(ZG, ViZG))         # Var(u | y), beta known
  C22 <- Cc + crossprod(XViZG, C11 %*% XViZG)
  V <- rbind(cbind(C11, C12), cbind(t(C12), C22)) * s2
  beta <- nlme::fixef(model)
  u <- as.vector(Matrix::crossprod(ViZG, R$y - X %*% beta))
  list(beta = beta, u = u, V = V, cond = diag(Cc) * s2, tab = do.call(rbind, unname(tab)), rows = rows, data = dat,
       fixed = pm_ranef_fixed_factors(stats::formula(model), dat), family = NULL,
       cond_label = "Var(u | y), fixed effects known (as lme4 condVar)")
}

## The factors of the fixed-effect formula (variables that are factors or
## character in the data).
pm_ranef_fixed_factors <- function(f, data) {
  v <- intersect(all.vars(stats::delete.response(stats::terms(f))), names(data))
  v[vapply(v, function(x) is.factor(data[[x]]) || is.character(data[[x]]), logical(1))]
}
