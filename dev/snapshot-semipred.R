# Snapshot semipred() output (result tables, plot layer data, plot pages) for
# option combinations written in the 1.1.1 spelling, rewritten to the 2.0
# spelling when the library has it (as vignette("migrating") says).
# Usage: Rscript --no-init-file dev/snapshot-semipred.R <library> <out.rds>
# then compare two snapshots with all.equal(a, b, tolerance = 0).
lib <- commandArgs(TRUE)[1]; out <- commandArgs(TRUE)[2]
suppressPackageStartupMessages({library(predictmeans, lib.loc = lib); library(ggplot2)})
set.seed(1)
n <- 240
d <- data.frame(x = runif(n, 0, 3), z = runif(n, 0, 2), f = factor(rep(c("a", "b"), each = n/2)))
d$y <- sin(2*d$x) + (d$f == "b")*d$x/2 + d$z + rnorm(n, sd = 0.3)
d$cnt <- rpois(n, exp(0.3 + sin(d$x)/2))
# the fits, in the 1.1.1 spelling (prt), as 2.0 when the library is 2.0
new_lib <- !"prt" %in% names(formals(semireg))
fit_call <- function(cl) {
  a <- as.list(cl)
  if (new_lib) a$prt <- NULL
  suppressMessages(suppressWarnings(eval(as.call(a))))
}
fits <- list(
  fit = quote(semireg(y ~ s(x, by = f) + s(z), data = d)),
  fit1 = quote(semireg(y ~ s(x) + z, data = d)),
  fp = quote(semireg(cnt ~ f + s(x), data = d, family = poisson())),
  f2 = quote(semireg(y ~ s(x, z, bs = "ztps", k = 6), data = d)),
  ft = quote(semireg_tmb(cnt ~ s(x, by = f, bs = "os", k = 6), family = glmmTMB::nbinom2, ziformula = ~1, data = d)),
  ftg = quote(semireg_tmb(y ~ s(x) + z, data = d)))
for (nm in names(fits)) assign(nm, fit_call(fits[[nm]]))
fit_summary <- function(m) list(logLik = as.numeric(logLik(m$semer)), fixef = unlist(lme4::fixef(m$semer)),
                                CovMat = as.matrix(m$CovMat), u = unlist(m$u_lst), class = class(m$semer))
xs <- seq(0.2, 2.8, length.out = 7)
calls <- list(
  one = quote(semipred(fit1, covariate = "x")),
  grp = quote(semipred(fit, "f", "x", covariateV = xs)),
  lev10 = quote(semipred(fit, "f", "x", covariateV = xs, level = 0.1)),
  sim = quote(semipred(fit, "f", "x", covariateV = xs, simultaneous = TRUE)),
  band = quote(semipred(fit, "f", "x", covariateV = xs, band = "simultaneous")),
  contr = quote(semipred(fit, "f", "x", covariateV = xs, contr = c("b", "a"))),
  pois = quote(semipred(fp, "f", "x", covariateV = xs)),
  pois_contr = quote(semipred(fp, "f", "x", covariateV = xs, contr = c(2, 1), contr_scale = "response")),
  deriv = quote(semipred(fit1, covariate = "x", covariateV = xs, deriv = 1)),
  plotset = quote(semipred(fit, "f", "x", trellis = FALSE, ci = FALSE, point = FALSE, jitterv = 0.1, prt = TRUE)),
  scales = quote(semipred(fit, "f", "x", scales = "free_y", plotord = 1)),
  surf = quote(semipred(f2, covariate = c("x", "z"))),
  named = quote(semipred(semireg = fit1, covariate = "x", covariateV = xs)),
  tmb = quote(semipred(ft, "f", "x", covariateV = xs, contr = c("b", "a"))),
  tmb_g = quote(semipred(ftg, covariate = "x", covariateV = xs)))
pm2 <- function(cl) {
  if (!"conf.level" %in% names(formals(semipred))) return(cl)
  a <- as.list(cl); n <- names(a)
  if ("level" %in% n) { a$conf.level <- 1 - a$level; a$level <- NULL }
  if ("covariateV" %in% n) { a$at <- a$covariateV; a$covariateV <- NULL }
  if ("semireg" %in% n) { a$model <- a$semireg; a$semireg <- NULL }
  if ("simultaneous" %in% n) { if (isTRUE(a$simultaneous)) a$band <- "simultaneous"; a$simultaneous <- NULL }
  if ("prt" %in% n) { a$plot <- a$prt; a$prt <- NULL }
  pl <- intersect(n, c("trellis", "scales", "plotord", "ci", "point", "jitterv", "threeD"))
  if (length(pl)) { a$plot_args <- as.call(c(as.name("list"), a[pl])); a[pl] <- NULL }
  as.call(a)
}
plotdata <- function(p) if (inherits(p, "ggplot")) lapply(seq_along(p$layers), function(i) ggplot2::layer_data(p, i)) else p
R <- list()
for (nm in names(calls)) {
  set.seed(42); f <- tempfile(fileext = ".pdf"); grDevices::pdf(f)
  r <- tryCatch(suppressMessages(suppressWarnings(eval(pm2(calls[[nm]])))),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  grDevices::dev.off()
  pages <- length(gregexpr("/Type /Page[^s]", paste(readLines(f, warn = FALSE), collapse = "\n"))[[1]])
  if (is.list(r)) {
    r <- unclass(r)
    r$settings$version <- NULL
    set.seed(44); r$plt <- plotdata(r$plt)
  }
  R[[nm]] <- list(result = r, pages = pages)
}
R$fits <- lapply(names(fits), function(nm) fit_summary(get(nm)))
saveRDS(R, out)
cat(length(R), "calls;", sum(vapply(R, function(x) is.character(x$result), logical(1))), "errors\n")
