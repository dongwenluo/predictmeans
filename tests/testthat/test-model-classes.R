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

test_that("a new class's coefficients are matched to the design by name, or the call stops", {
  # an empty cell aliases A2:B2; the made-up class lists its coefficients in
  # reverse order, the aliased one as NA (as lm() does)
  d <- expand.grid(A = factor(1:3), B = factor(1:2), r = 1:4)
  d <- d[!(d$A == "2" & d$B == "2"), ]
  d$y <- c(10, 12, 15)[d$A] + 2 * (d$B == "2") + sin(seq_len(nrow(d)))
  fl <- lm(y ~ A * B, d)
  rev_fit <- structure(list(fit = fl, call = quote(lm(y ~ A * B, d))), class = "pm_revmodel")
  registerS3method("coef", "pm_revmodel", function(object, ...) rev(coef(object$fit)))
  registerS3method("vcov", "pm_revmodel", function(object, ...) { v <- vcov(object$fit); v[rev(rownames(v)), rev(colnames(v))] })
  registerS3method("terms", "pm_revmodel", function(x, ...) terms(x$fit))
  registerS3method("df.residual", "pm_revmodel", function(object, ...) object$fit$df.residual)
  registerS3method("model.frame", "pm_revmodel", function(formula, ...) model.frame(formula$fit))
  q <- function(e) suppressWarnings(suppressMessages(e))
  a <- q(predictmeans(rev_fit, "A:B", pairwise = TRUE, data = d))
  b <- q(predictmeans(fl, "A:B", pairwise = TRUE))
  expect_equal(a$means$mean, b$means$mean, tolerance = 1e-10)
  expect_equal(a$means$std.error, b$means$std.error, tolerance = 1e-10)
  expect_equal(a$pairs$estimate, b$pairs$estimate, tolerance = 1e-10)
  # coefficient names that are not columns of the design: stop, never guess
  bad_fit <- structure(list(fit = fl, call = quote(lm(y ~ A * B, d))), class = "pm_badmodel")
  registerS3method("coef", "pm_badmodel", function(object, ...) { b <- coef(object$fit); names(b) <- paste0("b", seq_along(b)); b })
  registerS3method("vcov", "pm_badmodel", function(object, ...) { v <- vcov(object$fit); dimnames(v) <- list(paste0("b", 1:5), paste0("b", 1:5)); v })
  registerS3method("terms", "pm_badmodel", function(x, ...) terms(x$fit))
  registerS3method("df.residual", "pm_badmodel", function(object, ...) object$fit$df.residual)
  registerS3method("model.frame", "pm_badmodel", function(formula, ...) model.frame(formula$fit))
  expect_error(q(predictmeans(bad_fit, "A", data = d)), "b1, b2, b3, b4, b5, ... have no column in the design", fixed = TRUE)
})
