## Level and power of the tests of variance components (varcomp(test = TRUE),
## R/vc-test.R): the default (method = "auto": the boundary mixtures, the exact
## restricted LRT for the only random effect of a gaussian fit) vs the
## parametric bootstrap (199 draws per test). Each scenario with the tested
## component 0 (level) and > 0 (power):
##   A gaussian, crossed (1 | a) + (1 | b): b's variance       (8 x 6 x 2)
##   B gaussian, (x | s): the slope variance                    (15 subjects x 6)
##   C gaussian, one-way with 4 levels (flagged): exact test    (4 x 6)
##   D binary, (1 | g), 20 clusters of 4 (flagged: Laplace)
##   E Poisson, (1 | g), 15 clusters of 6
##   F negative binomial (glmmTMB nbinom2), (1 | g), 15 x 6
##   G semireg s(x, f1*f2): the f2 deviation curves
##   Rscript dev/sim-vc-defaults.R [scenarios, e.g. ABC] [nrep] [outdir]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
which_sc <- if (length(args)) strsplit(args[1], "")[[1]] else LETTERS[1:7]
nrep <- if (length(args) > 1) as.integer(args[2]) else 200
outdir <- if (length(args) > 2) args[3] else tempdir()
q <- function(e) suppressWarnings(suppressMessages(e))
pv <- function(fit, row) {
  a <- q(pm_vc_test(fit, rows = row, nsim = 4000))
  b <- q(pm_vc_test(fit, rows = row, method = "bootstrap", nboot = 199))
  c(default = a$p.value[1], bootstrap = b$p.value[1])
}
sc <- list(
  A = function(seed, alt) {
    set.seed(seed)
    d <- expand.grid(a = factor(1:8), b = factor(1:6), r = 1:2)
    d$y <- rnorm(8, 0, 0.7)[d$a] + rnorm(6, 0, if (alt) 0.5 else 0)[d$b] + rnorm(nrow(d))
    pv(q(lme4::lmer(y ~ 1 + (1 | a) + (1 | b), d)), "b.(Intercept)")
  },
  B = function(seed, alt) {
    set.seed(seed)
    d <- expand.grid(x = 0:5, s = factor(1:15))
    d$y <- 1 + 0.5*d$x + rnorm(15, 0, 1)[d$s] + rnorm(15, 0, if (alt) 0.3 else 0)[d$s]*d$x + rnorm(nrow(d), 0, 0.7)
    pv(q(lme4::lmer(y ~ x + (x | s), d)), "s.x")
  },
  C = function(seed, alt) {
    set.seed(seed)
    d <- expand.grid(g = factor(1:4), r = 1:6)
    d$y <- rnorm(4, 0, if (alt) 0.7 else 0)[d$g] + rnorm(nrow(d))
    pv(q(lme4::lmer(y ~ 1 + (1 | g), d)), "g.(Intercept)")
  },
  D = function(seed, alt) {
    set.seed(seed)
    d <- expand.grid(r = 1:4, g = factor(1:20))
    d$y <- rbinom(nrow(d), 1, plogis(-0.3 + rnorm(20, 0, if (alt) 0.8 else 0)[d$g]))
    pv(q(lme4::glmer(y ~ 1 + (1 | g), family = binomial, data = d)), "g.(Intercept)")
  },
  E = function(seed, alt) {
    set.seed(seed)
    d <- expand.grid(r = 1:6, g = factor(1:15))
    d$y <- rpois(nrow(d), exp(1 + rnorm(15, 0, if (alt) 0.4 else 0)[d$g]))
    pv(q(lme4::glmer(y ~ 1 + (1 | g), family = poisson, data = d)), "g.(Intercept)")
  },
  F = function(seed, alt) {
    set.seed(seed)
    d <- expand.grid(r = 1:6, g = factor(1:15))
    d$y <- rnbinom(nrow(d), size = 3, mu = exp(1 + rnorm(15, 0, if (alt) 0.4 else 0)[d$g]))
    pv(q(glmmTMB::glmmTMB(y ~ 1 + (1 | g), family = glmmTMB::nbinom2, data = d)), "g.(Intercept)")
  },
  G = function(seed, alt) {
    set.seed(seed)
    d <- expand.grid(x = seq(0, 1, length = 25), f1 = factor(c("p", "q")), f2 = factor(c("u", "v", "w")), rep = 1:4)
    d$id <- interaction(d$f1, d$f2, d$rep)
    d$y <- sin(2*pi*d$x) + 1.5*(d$f1 == "q")*sin(pi*d$x)^2 + if (alt) 1.0*(d$f2 == "w")*sin(pi*d$x)^2 else 0
    d$y <- d$y + rnorm(nlevels(d$id), 0, 0.2)[d$id] + rnorm(nrow(d), 0, 0.25)
    pv(q(semireg(y ~ s(x, f1*f2, bs = "os", k = 8) + (1 | id), data = d, prt = FALSE)), "s_x_f2.(Intercept)")
  })
for (s in which_sc) {
  logf <- file.path(outdir, paste0("sim-vc-", s, ".log"))
  t0 <- proc.time()[["elapsed"]]
  res <- lapply(c(null = FALSE, alt = TRUE), function(alt) {
    n <- if (alt) ceiling(nrep/2) else nrep
    r <- parallel::mclapply(seq_len(n), function(i) {
      v <- tryCatch(sc[[s]](10000*match(s, LETTERS) + i + if (alt) 5000 else 0, alt), error = function(e) c(NA, NA))
      cat(s, alt, i, v, "\n", file = logf, append = TRUE)
      v
    }, mc.cores = parallel::detectCores())
    do.call(rbind, r)
  })
  saveRDS(res, file.path(outdir, paste0("sim-vc-", s, ".rds")))
  f <- function(m) { m <- m[stats::complete.cases(m), , drop = FALSE]; sprintf("default %.3f, bootstrap %.3f (n %d)",
                                                                              mean(m[, 1] < 0.05), mean(m[, 2] < 0.05), nrow(m)) }
  cat(sprintf("%s: level %s | power %s | %.0f s\n", s, f(res$null), f(res$alt), proc.time()[["elapsed"]] - t0))
}
