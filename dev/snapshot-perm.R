# Snapshot of the permutation and bootstrap functions (results and what they
# print), with calls written in the 1.1.1 spelling and rewritten to the 2.0
# spelling when the library has it (as vignette("migrating") says).
# Usage: Rscript --no-init-file dev/snapshot-perm.R <library> <out.rds>
# then compare two snapshots with all.equal(a, b, tolerance = 0).
lib <- commandArgs(TRUE)[1]; out <- commandArgs(TRUE)[2]
suppressPackageStartupMessages(library(predictmeans, lib.loc = lib))
Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
mt <- mtcars; mt$cyl <- factor(mt$cyl)
m_lm <- lm(mpg ~ wt + cyl, data = mt)
m_lmer <- lme4::lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = Oats)
m_lmer0 <- lme4::lmer(yield ~ nitro + Variety + (1 | Block), data = Oats)
data(cbpp, package = "lme4")
m_glmer <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = cbpp)
m_glmer0 <- lme4::glmer(cbind(incidence, size - incidence) ~ 1 + (1 | herd), family = binomial, data = cbpp)
m_glm <- glm(cbind(incidence, size - incidence) ~ period, family = binomial, data = cbpp)
f_obj <- mpg ~ wt + cyl
m_fo <- lm(f_obj, data = mt)
new <- !"prt" %in% names(formals(permmodels))
calls <- list(
  lm_fl = quote(permmodels(m_lm, nperm = 99, type = "II", seed = 1, prt = FALSE)),
  lm_raw = quote(permmodels(m_lm, nperm = 99, seed = 1, method = "raw", prt = FALSE)),
  lm_robust = quote(permmodels(m_lm, nperm = 99, seed = 1, robust = TRUE, prt = FALSE)),
  lmer_fl = quote(permmodels(m_lmer, nperm = 19, seed = 2, prt = FALSE)),
  lmer_raw = quote(permmodels(m_lmer, nperm = 19, seed = 2, method = "raw", prt = FALSE)),
  glm_fl = quote(permmodels(m_glm, nperm = 49, seed = 3, prt = FALSE)),
  fo = quote(permmodels(m_fo, nperm = 49, seed = 4, method = "raw", fo = f_obj, prt = FALSE)),
  boot = quote(bootmodels(m_lmer, nsim = 19, seed = 5, prt = FALSE)),
  permlmer = quote(permlmer(m_lmer0, m_lmer, nperm = 19, seed = 6)),
  permlmer_named = quote(permlmer(lmer0 = m_lmer0, lmer1 = m_lmer, nperm = 19, seed = 6)),
  permglmer = quote(permglmer(m_glmer0, m_glmer, nperm = 19, seed = 7)),
  permstats = quote(permstats(permmodels(m_lm, nperm = 99, type = "II", seed = 1, prt = FALSE), level = 0.1)))
pm2 <- function(cl) {
  if (!new) return(cl)
  a <- as.list(cl); n <- names(a)
  a$prt <- NULL; a$fo <- NULL
  if (identical(a[[1]], as.name("permglmer"))) a[[1]] <- as.name("permlmer")
  if ("lmer0" %in% n) { a$model0 <- a$lmer0; a$lmer0 <- NULL }
  if ("lmer1" %in% n) { a$model1 <- a$lmer1; a$lmer1 <- NULL }
  if ("level" %in% n) { a$conf.level <- 1 - a$level; a$level <- NULL }
  if (identical(a[[1]], as.name("permstats"))) a[[2]] <- pm2(a[[2]])
  as.call(a)
}
R <- list()
for (nm in names(calls)) {
  printed <- NULL
  r <- tryCatch({
    res <- suppressWarnings(suppressMessages(eval(pm2(calls[[nm]]))))
    printed <- capture.output(print(res))
    res
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.list(r)) { r <- unclass(r); r$notes <- NULL }
  R[[nm]] <- list(result = r, printed = printed)
}
saveRDS(R, out)
cat(length(R), "calls;", sum(vapply(R, function(x) is.character(x$result), logical(1))), "errors\n")
