## Permutation tests of nested GLMMs (glmer or glmmTMB), called by
## permlmer() and its aliases permglmer() and permglmmTMB().
## dev/permutation-plan.md, sections 4 and 7; dev/sim-study/ for the
## simulations behind the choices.
##
## "working" (route A): the H0 model is linearised (pm_working()); its
##   working residuals are decorrelated with the order-invariant square root
##   of V0 = A A' + diag(resvar), permuted, correlated again and added to
##   X beta0: z* = X b0 + W e[pi]. The statistic is the likelihood-ratio
##   statistic of weighted linear mixed models (lmer, weights 1/resvar)
##   fitted to z* (REML when the fixed effects are the same, i.e. a test of
##   random effects; ML otherwise), and the observed one the same from z.
##   For random effects this is the Lee-Braun test on the working scale, for
##   fixed effects Freedman-Lane (H0 fitted values kept). In simulations it
##   kept the level for counts and binomial data with several trials.
## "quantile" (route B): randomised quantile residuals of the 0/1 responses
##   given the H0 conditional means are permuted and mapped back to 0/1 with
##   the same means; both GLMMs are refitted (glmer with nAGQ = 0, observed
##   and permuted alike). The working scale fails for 0/1 responses with rare
##   events or tiny clusters (36% at the 5% level); this kept it (5.6-6.2%).
pm_permglmm <- function(m0, m1, nperm, ncore, seed, perms, method, plot, labels) {
  tmb <- inherits(m0, "glmmTMB")
  if (tmb != inherits(m1, "glmmTMB") || (!tmb && !(inherits(m0, "glmerMod") && inherits(m1, "glmerMod")))) {
    stop("Both models must be glmer models, or both glmmTMB models.", call. = FALSE)
  }
  W0 <- pm_working(m0)
  W1y <- pm_working(m1)$y
  if (length(W0$y) != length(W1y) || any(W0$y != W1y)) stop("Please check the response in your model!", call. = FALSE)
  f0 <- stats::formula(m0)
  f1 <- stats::formula(m1)
  fe0 <- colnames(W0$X)
  fe1 <- colnames(pm_working(m1)$X)
  bars <- function(f) vapply(pm_findbars(f), function(b) paste(deparse(b), collapse = ""), "")
  if (!all(fe0 %in% fe1) || !all(bars(f0) %in% bars(f1))) {
    stop("The model ", labels[1], " must be nested within the model ", labels[2], ".", call. = FALSE)
  }
  ref <- setequal(fe0, fe1)              # same fixed effects: a test of random effects
  if (method == "auto") method <- if (W0$binary) "quantile" else "working"
  if (method == "quantile" && !W0$binary) {
    stop("method = \"quantile\" is for 0/1 responses; use method = \"working\".", call. = FALSE)
  }
  n <- length(W0$z)
  if (!is.null(perms)) nperm <- ncol(as.matrix(perms))
  idx <- pm_perm_source(n, seed, perms)(nperm)
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  lrt <- function(a1, a0) as.numeric(2 * (stats::logLik(a1) - stats::logLik(a0)))

  if (method == "working") {
    olre <- vapply(pm_re_terms(m1), function(tm) tm$nl >= n, logical(1))
    if (any(olre)) {
      stop("A random effect with one level per observation cannot be used on the working scale (it is not separable ",
           "from the working residual); a parametric bootstrap (bootmodels()) can test it.", call. = FALSE)
    }
    d <- pm_model_rows(m1)              # the larger model's data has every variable
    d$.pm_z <- W0$z
    ctl <- lme4::lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")
    fit_w <- function(f) {
      # the weights are found in the formula's environment, this frame
      .pm_wts <- 1 / W0$resvar
      wf <- pm_working_formula(f)
      environment(wf) <- environment()
      out <- tryCatch(quiet(lme4::lmer(wf, data = d, weights = .pm_wts, REML = ref, control = ctl)),
                      error = function(e) e)
      if (inherits(out, "error")) {
        stop("The working linear mixed model could not be fitted by lmer (", conditionMessage(out),
             "); random effects must be terms that lmer can fit.", call. = FALSE)
      }
      out
    }
    l0 <- fit_w(f0)
    l1 <- fit_w(f1)
    W <- pm_whiten(W0$A, W0$resvar)
    e <- W$solve(W0$z - W0$xb)
    obs <- lrt(l1, l0)
    stat <- function(k) {
      zs <- W0$xb + W$mult(e[idx[, k]])
      lrt(quiet(lme4::refit(l1, zs)), quiet(lme4::refit(l0, zs)))
    }
  } else {
    # 0/1 responses: refits of the GLMMs themselves
    if (tmb) {
      d <- pm_model_rows(m1)
      yname <- all.vars(f0[[2]])
      if (length(yname) != 1) stop("method = \"quantile\" needs a response that is one variable.", call. = FALSE)
      g0 <- m0
      g1 <- m1
      refit_y <- function(m, y) {
        dd <- d
        dd[[yname]] <- if (is.factor(d[[yname]])) factor(levels(d[[yname]])[y + 1], levels = levels(d[[yname]])) else y
        quiet(pm_refit(m, data = dd))
      }
    } else {
      ctl <- lme4::glmerControl(optimizer = "bobyqa", calc.derivs = FALSE, check.conv.singular = "ignore")
      g0 <- quiet(pm_refit(m0, nAGQ = 0, control = ctl))
      g1 <- quiet(pm_refit(m1, nAGQ = 0, control = ctl))
      refit_y <- function(m, y) quiet(lme4::refit(m, y))
    }
    mu0 <- if (tmb) stats::fitted(g0) else g0@resp$mu      # conditional on the H0 random effects
    y <- W0$y
    lo <- stats::pbinom(y - 1, 1, mu0)
    hi <- stats::pbinom(y, 1, mu0)
    u <- if (is.null(seed)) stats::runif(n) else pm_with_seed(seed + 1, stats::runif(n))
    q <- stats::qnorm(pmin(pmax(lo + u * (hi - lo), 1e-12), 1 - 1e-12))
    obs <- lrt(g1, g0)
    stat <- function(k) {
      ys <- stats::qbinom(stats::pnorm(q[idx[, k]]), 1, mu0)
      lrt(refit_y(g1, ys), refit_y(g0, ys))
    }
  }

  sims <- pm_perm_drop(pm_perm_map(seq_len(nperm), stat, ncore, packages = c("lme4", "glmmTMB")))
  s <- unlist(sims$res)
  ok <- !is.na(s)
  failures <- sims$failures + sum(!ok)
  s <- round(pmax(s[ok], 0), 4)
  obs <- round(max(obs, 0), 4)
  aod <- quiet(stats::anova(m0, m1))
  aod$`Perm-p` <- c(NA, pm_perm_p(obs, s, digits = 4))
  attr(aod, "method") <- method
  if (failures > 0) attr(aod, "failures") <- failures
  if (plot) {
    grDevices::dev.new()
    graphics::plot(stats::density(c(obs, s), kernel = "epanechnikov"), col = "blue", lwd = 2, xlab = "",
                   main = "LR Test's density kernels")
    graphics::abline(v = obs, col = "red")
  }
  aod
}
