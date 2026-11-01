# Snapshot predictmeans() output for ~30 option combinations (result object,
# ggplot layer data, printed output, graphics pages), to check that a refactor
# changes nothing. Usage:
#   Rscript --no-init-file dev/snapshot-predictmeans.R <library> <out.rds>
# then compare two snapshots with all.equal(a, b, tolerance = 0).
lib <- commandArgs(TRUE)[1]; out <- commandArgs(TRUE)[2]
suppressPackageStartupMessages({library(predictmeans, lib.loc = lib); library(ggplot2)})
Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
Om <- Oats[!(Oats$nitro == "0.6" & Oats$Variety == "Victory"), ]
data(Drug); data(Clinical); data(ATP)
m_lm <- lm(yield ~ Block + nitro*Variety, Oats)
m_lmer <- lmerTest::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
m_lme <- lme(yield ~ nitro*Variety, random = ~1|Block/Variety, data = Oats)
m_glm <- glm(cbind(R, N - R) ~ Drug + log10Dose, binomial, Drug)
m_glmer <- glmer(cbind(Favorable, Total - Favorable) ~ Treatment + (1|Clinic), family = binomial, data = Clinical)
m_miss <- lmerTest::lmer(yield ~ nitro*Variety + (1|Block), data = Om)
m_atp <- lmerTest::lmer(ATP ~ A*B*time + (1|heart), data = ATP)
m_log <- lm(log(yield) ~ nitro*Variety, Oats)
m_ml <- lmerTest::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats, REML = FALSE)
set.seed(3); dp <- expand.grid(A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")), rep = 1:5)
dp$y <- rpois(nrow(dp), exp(1 + as.numeric(dp$A)/2 + as.numeric(dp$B)/4))
m_pois <- glm(y ~ A*B, family = poisson, data = dp)
Oc <- Oats; levels(Oc$Variety) <- c("Golden:Rain", "Marv", "Vic")
m_colon <- lm(yield ~ Block + nitro*Variety, data = Oc)
perm <- suppressWarnings(permmodels(m_lmer, nperm = 19, prt = FALSE, ncore = 1, seed = 1))
calls <- list(
  lm_def = quote(predictmeans(m_lm, "nitro")),
  lm_int = quote(predictmeans(m_lm, "nitro:Variety")),
  lm_pw = quote(predictmeans(m_lm, "nitro:Variety", pairwise = TRUE)),
  lm_BH = quote(predictmeans(m_lm, "nitro:Variety", adj = "BH")),
  lm_tukey = quote(predictmeans(m_lm, "nitro", adj = "tukey")),
  lm_bonf = quote(predictmeans(m_lm, "nitro:Variety", adj = "bonferroni")),
  lm_atvar = quote(predictmeans(m_lm, "nitro:Variety", atvar = "Variety", adj = "BH")),
  lm_atvar_tk = quote(predictmeans(m_lm, "nitro:Variety", atvar = "Variety", adj = "tukey")),
  lm_letterCI = quote(predictmeans(m_lm, "nitro", letterCI = TRUE)),
  lm_decr = quote(predictmeans(m_lm, "nitro:Variety", adj = "BH", meandecr = TRUE)),
  lm_allplots = quote(predictmeans(m_lm, "nitro:Variety", which.plots = c("mean", "bar", "pval", "back"))),
  lm_plotord = quote(predictmeans(m_lm, "nitro:Variety", plotord = c(2, 1), lsd_bar = FALSE)),
  lm_bar1 = quote(predictmeans(m_lm, "nitro", which.plots = c("mean", "bar"), lineplot = FALSE)),
  lm_noplot = quote(predictmeans(m_lm, "nitro:Variety", adj = "BH", plot = FALSE)),
  log_trans = quote(predictmeans(m_log, "nitro", trans = exp, adj = "BH")),
  log_trans_off = quote(predictmeans(m_log, "nitro", trans = exp, transOff = 1)),
  lmer_int = quote(predictmeans(m_lmer, "nitro:Variety", adj = "BH")),
  lmer_unadj = quote(predictmeans(m_lmer, "nitro:Variety", adj = "tukey", vcov_adjust = FALSE)),
  lmer_atvar = quote(predictmeans(m_lmer, "nitro:Variety", atvar = "Variety", adj = "tukey")),
  lmer_perm = quote(predictmeans(m_lmer, "nitro", permlist = perm, adj = "BH", ncore = 1)),
  lmer_perm_at = quote(predictmeans(m_lmer, "nitro:Variety", permlist = perm, atvar = "Variety", ncore = 1)),
  lmer_ml = quote(suppressMessages(predictmeans(m_ml, "nitro", adj = "BH"))),
  lme_int = quote(predictmeans(m_lme, "nitro:Variety", adj = "holm")),
  glm_def = quote(predictmeans(m_glm, "Drug", pairwise = TRUE)),
  glm_cov = quote(predictmeans(m_glm, "Drug", covariate = 0.5, adj = "BH")),
  glmer = quote(predictmeans(m_glmer, "Treatment", pairwise = TRUE)),
  miss = quote(suppressWarnings(predictmeans(m_miss, "nitro:Variety", adj = "BH"))),
  miss_marg = quote(suppressWarnings(predictmeans(m_miss, "nitro"))),
  atp_3way = quote(predictmeans(m_atp, "A:B:time", atvar = "time", adj = "BH")),
  atp_2way = quote(predictmeans(m_atp, "A:time", plotord = c(2, 1))),
  user_df = quote(predictmeans(m_lmer, "nitro", Df = 20, adj = "BH")),
  weights = quote(predictmeans(m_lm, "nitro", weights = "proportional")),
  pois_count = quote(predictmeans(m_pois, "A", trans = exp, count = TRUE, pairwise = TRUE)),
  pois_atvar = quote(predictmeans(m_pois, "A:B", atvar = "B", adj = "BH", trans = exp)),
  pois_letterCI = quote(predictmeans(m_pois, "A:B", pairwise = TRUE, letterCI = TRUE, trans = exp, meandecr = FALSE)),
  colon_levels = quote(predictmeans(m_colon, "nitro:Variety", atvar = "nitro", adj = "BH")))
plotdata <- function(p) if (inherits(p, "ggplot")) lapply(seq_along(p$layers), function(i) ggplot2::layer_data(p, i)) else p
R <- list()
for (nm in names(calls)) {
  set.seed(42); f <- tempfile(fileext = ".pdf"); grDevices::pdf(f)
  printed <- NULL
  r <- tryCatch({
    res <- suppressWarnings(eval(calls[[nm]]))
    printed <- capture.output(print(res))
    res
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  grDevices::dev.off()
  pages <- length(gregexpr("/Type /Page[^s]", paste(readLines(f, warn = FALSE), collapse = "\n"))[[1]])
  if (is.list(r)) {
    r$settings$version <- NULL
    for (pn in c("predictmeansPlot", "predictmeansBarPlot")) {
      set.seed(44); if (!is.null(r[[pn]])) r[[pn]] <- if (inherits(r[[pn]], "ggplot")) plotdata(r[[pn]]) else lapply(r[[pn]], plotdata)
    }
  }
  set.seed(43); R[[nm]] <- list(result = unclass(r), printed = printed, pages = pages)
}
saveRDS(R, out)
cat(length(R), "calls;", sum(vapply(R, function(x) is.character(x$result), logical(1))), "errors\n")
