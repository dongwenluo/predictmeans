## Every main function on every model class predictmeans() handles. A
## function either works, or refuses the class with a message that says so
## (the refusals listed in 'refused'). Added 2026-11-17 so that a class gap
## cannot come back unnoticed.
test_that("every main function works on the 13 model classes, or refuses with a clear message", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("mgcv")
  set.seed(1)
  d <- expand.grid(A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")), Blk = factor(1:8), r = 1:2)
  d$x <- runif(nrow(d), 0, 2)
  u <- rnorm(8, sd = 0.5)
  d$y <- 2 + as.numeric(d$A)*0.5 + 0.8*d$x + 0.4*(d$A == "a3")*d$x + u[d$Blk] + rnorm(nrow(d), sd = 0.6)
  d$cnt <- rpois(nrow(d), exp(0.2 + 0.3*as.numeric(d$A) + 0.4*d$x + u[d$Blk]))
  d$bin <- rbinom(nrow(d), 1, plogis(-1 + 0.5*as.numeric(d$A) + 0.6*d$x + u[d$Blk]))
  q <- function(e) suppressWarnings(suppressMessages(e))
  mods <- q(list(
    lm = lm(y ~ A*x + B, d), aov = aov(y ~ A*x + B, d), aovlist = aov(y ~ A*x + B + Error(Blk), d),
    glm_pois = glm(cnt ~ A*x + B, poisson, d), glm_bin = glm(bin ~ A*x + B, binomial, d),
    gls = nlme::gls(y ~ A*x + B, correlation = nlme::corCompSymm(form = ~1 | Blk), data = d),
    lme = nlme::lme(y ~ A*x + B, random = ~1 | Blk, data = d),
    lmer = lme4::lmer(y ~ A*x + B + (1 | Blk), d),
    glmer_pois = lme4::glmer(cnt ~ A*x + B + (1 | Blk), d, family = poisson),
    glmer_bin = lme4::glmer(bin ~ A*x + B + (1 | Blk), d, family = binomial),
    tmb_pois = glmmTMB::glmmTMB(cnt ~ A*x + B + (1 | Blk), d, family = poisson),
    tmb_nb = glmmTMB::glmmTMB(cnt ~ A*x + B + (1 | Blk), d, family = glmmTMB::nbinom2),
    gam = mgcv::gam(y ~ A*x + B + s(Blk, bs = "re"), data = d, method = "REML")))
  cm <- rbind(c(1, -1, 0), c(0, 1, -1))
  fns <- list(
    predictmeans = function(m) predictmeans(m, "A", pairwise = TRUE),
    covariatemeans = function(m) covariatemeans(m, "A", covariate = "x"),
    cov_contr = function(m) covariatemeans(m, "A", covariate = "x", contr = c("a1", "a3")),
    contrastmeans = function(m) contrastmeans(m, "A", ctrmatrix = cm),
    pm_joint = function(m) pm_joint(m, "A"),
    pm_slopes = function(m) pm_slopes(m, "x", "A"),
    varcomp = function(m) varcomp(m),
    R2_model = function(m) R2_model(m),
    residplot = function(m) residplot(m),
    pm_check = function(m) pm_check(m),
    CookD = function(m) CookD(m, group = "Blk", plot = FALSE),
    permmodels = function(m) permmodels(m, nperm = 5, ncore = 1, seed = 1),
    bootmodels = function(m) bootmodels(m, nsim = 5, ncore = 1, seed = 1),
    autoplot = function(m) ggplot2::autoplot(predictmeans(m, "A", pairwise = TRUE),
                                             interval = "comparison"),
    Kmatrix = function(m) Kmatrix(m, "A", list(x = c(0, 1))))
  refused <- c("CookD:aovlist", "permmodels:aovlist", "bootmodels:aovlist", "bootmodels:gam")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (f in names(fns)) for (nm in names(mods)) {
    r <- tryCatch(q(utils::capture.output(fns[[f]](mods[[nm]]))), error = function(e) e)
    key <- paste(f, nm, sep = ":")
    if (key %in% refused) {
      expect_s3_class(r, "error")
      expect_match(conditionMessage(r), "does not handle", label = key)
    } else {
      expect_false(inherits(r, "error"), label = paste(key, if (inherits(r, "error")) conditionMessage(r)))
    }
  }
})
