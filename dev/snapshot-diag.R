# Snapshot of the plot, diagnostic and small functions (results; plots as
# ggplot layer data), with calls written in the 1.1.1 spelling and rewritten
# to the 2.0 spelling when the library has it (vignette("migrating")).
# Usage: Rscript --no-init-file dev/snapshot-diag.R <library> <out.rds>
lib <- commandArgs(TRUE)[1]; out <- commandArgs(TRUE)[2]
suppressPackageStartupMessages({library(predictmeans, lib.loc = lib); library(ggplot2)})
Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
m_lmer <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
m_lm <- lm(yield ~ Block + nitro * Variety, data = Oats)
data(Clinical)
m_glmer <- lme4::glmer(cbind(Favorable, Total - Favorable) ~ Treatment + (1 | Clinic), family = binomial, data = Clinical)
pm <- predictmeans(m_lm, "nitro:Variety", adj = "BH")
new <- "re_level" %in% names(formals(residplot))
calls <- list(
  resid = quote(residplot(m_lmer, group = "nitro", plot = FALSE)),
  resid_l2 = quote(residplot(m_lmer, level = 2, plot = FALSE)),
  cook = quote(CookD(m_lm, plot = FALSE, idn = 5)),
  cook_g = quote(CookD(m_lmer, group = "Block", plot = FALSE)),
  pmplot = quote(PMplot(pm$p_valueMatrix, level = 0.1)),
  pmplot_m = quote(PMplot(pm$p_valueMatrix, mtitle = "T", margin = 3, legendx = 0.8)),
  vc = quote(varcomp(m_lmer, level = 0.9)),
  vc_glmer = quote(varcomp(m_glmer)),
  r2_lmer = quote(R2_model(m_lmer)),
  r2_glmer = quote(R2_model(m_glmer)),
  r2glmm_lmer = quote(R2_glmm(m_lmer)),
  r2glmm_glmer = quote(R2_glmm(m_glmer)))
pm2 <- function(cl) {
  if (!new) return(cl)
  a <- as.list(cl); n <- names(a); f <- as.character(a[[1]])
  if (f == "residplot" && "level" %in% n) { a$re_level <- a$level; a$level <- NULL }
  if (f == "CookD" && "idn" %in% n) { a$id <- a$idn; a$idn <- NULL }
  if (f == "PMplot") { if ("level" %in% n) { a$conf.level <- 1 - a$level; a$level <- NULL }; a$margin <- NULL; a$legendx <- NULL }
  if (f == "varcomp" && "level" %in% n) { a$conf.level <- a$level; a$level <- NULL }
  if (f == "R2_glmm") a[[1]] <- as.name("R2_model")
  as.call(a)
}
plotdata <- function(p) if (inherits(p, "ggplot")) lapply(seq_along(p$layers), function(i) ggplot2::layer_data(p, i)) else p
R <- list()
for (nm in names(calls)) {
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f)
  r <- tryCatch(suppressWarnings(suppressMessages(eval(pm2(calls[[nm]])))), error = function(e) paste("ERROR:", conditionMessage(e)))
  grDevices::dev.off()
  if (inherits(r, "ggplot")) r <- plotdata(r)
  if (is.list(r) && !is.data.frame(r)) r <- lapply(unclass(r), function(p) if (inherits(p, "ggplot")) plotdata(p) else p)
  R[[nm]] <- r
}
saveRDS(R, out)
cat(length(R), "calls;", sum(vapply(R, function(x) is.character(x) && length(x) == 1 && grepl("^ERROR", x), logical(1))), "errors\n")
