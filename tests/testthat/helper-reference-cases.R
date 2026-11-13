# Calls whose results are frozen in _data/reference-v1.1.1.rds. The reference
# was produced by running reference_cases() with predictmeans 1.1.1 (see
# dev/make-reference.R), so these tests check that later versions reproduce
# the 1.1.1 numbers for everything 1.1.1 supported.

reference_cases <- function() {
  R <- list()
  cap <- function(nm, expr) {
    R[[nm]] <<- tryCatch(suppressMessages(suppressWarnings(eval.parent(substitute(expr)))),
                         error = function(e) paste("ERROR:", conditionMessage(e)))
  }
  km <- function(...) {
    k <- Kmatrix(...)
    list(K = unclass(k$K)[, , drop = FALSE], fctnames = k$fctnames, response = k$response)
  }
  # 1.1.1 used unadjusted SEs, nlme's df for lme/gls and the residual df for
  # GLMMs; newer versions need vcov_adjust = FALSE and df_method =
  # "containment" (lme, gls) or "residual" (GLMMs) for that
  old_df <- function(f, m) {
    a <- names(formals(f))
    c(if ("vcov_adjust" %in% a) list(vcov_adjust = FALSE),
      if ("df_method" %in% a && inherits(m, c("lme", "gls"))) list(df_method = "containment"),
      if ("df_method" %in% a && inherits(m, c("glmerMod", "glmmTMB"))) list(df_method = "residual"))
  }
  cvm <- function(...) do.call(covariatemeans, c(list(...), old_df(covariatemeans, list(...)[[1]])))
  pmo <- function(...) {
    extra <- old_df(predictmeans, list(...)[[1]])
    p <- do.call(predictmeans, c(list(...), plot = FALSE, prtnum = FALSE, extra))
    # legacy slots only (the stable means/pairs/settings slots are newer)
    p[!grepl("Plot|plot", names(p)) & !names(p) %in% c("means", "pairs", "settings")]
  }

  Oats <- as.data.frame(nlme::Oats)
  fm_lme_num <- lme(yield ~ nitro*Variety, random = ~1|Block/Variety, data = Oats)
  cap("ex_K1", km(fm_lme_num, "Variety"))
  cap("ex_K2", km(fm_lme_num, "Variety", 0.5))
  cap("ex_K3", km(fm_lme_num, "Variety", "nitro", covariateV = seq(0, 0.6, 0.1)))
  cap("ex_K4", km(fm_lme_num, "Variety", "nitro"))
  cap("ex_cov1", cvm(fm_lme_num, "Variety", covariate = "nitro")$pltdf)
  cap("ex_cov2", cvm(fm_lme_num, "Variety", covariate = "nitro", covariateV = seq(0, 0.6, 0.1))$pltdf)

  Oats$nitro <- factor(Oats$nitro); Oats$Whole_plot <- Oats$Variety
  fm_lmer <- lmerTest::lmer(yield ~ nitro*Variety + (1|Block/Whole_plot), data = Oats)
  fm_lme <- lme(yield ~ nitro*Variety, random = ~1|Block/Variety, data = Oats)
  fm_lm <- lm(yield ~ Block + nitro*Variety, data = Oats)
  fm_aov <- aov(yield ~ nitro*Variety + Error(Block/Variety), data = Oats)
  for (t in c("nitro", "Variety", "nitro:Variety")) {
    cap(paste0("oats_lmer_K_", t), km(fm_lmer, t))
    cap(paste0("oats_lme_K_", t), km(fm_lme, t))
    cap(paste0("oats_lm_K_", t), km(fm_lm, t))
  }
  cap("oats_pm_lmer_nitro", pmo(fm_lmer, "nitro", adj = "BH"))
  cap("oats_pm_lmer_int", pmo(fm_lmer, "nitro:Variety", adj = "BH"))
  cap("oats_pm_lmer_atvar", pmo(fm_lmer, "nitro:Variety", atvar = "Variety", adj = "BH"))
  cap("oats_pm_lme_int", pmo(fm_lme, "nitro:Variety", pairwise = TRUE))
  cap("oats_pm_lm_int", pmo(fm_lm, "nitro:Variety", adj = "tukey"))
  cap("oats_pm_aov", pmo(fm_aov, "nitro"))

  Om <- Oats[!(Oats$nitro == "0.6" & Oats$Variety == "Victory"), ]
  fm_miss <- lm(yield ~ nitro*Variety, data = Om)
  fm_miss_lmer <- lmerTest::lmer(yield ~ nitro*Variety + (1|Block), data = Om)
  cap("miss_lm_K_int", km(fm_miss, "nitro:Variety"))
  cap("miss_lm_K_nitro", km(fm_miss, "nitro"))
  cap("miss_lm_K_var", km(fm_miss, "Variety"))
  cap("miss_lmer_K_int", km(fm_miss_lmer, "nitro:Variety"))
  cap("miss_pm", pmo(fm_miss, "nitro:Variety", pairwise = TRUE))

  data("ATP", package = "predictmeans", envir = environment())
  fm_atp <- lmerTest::lmer(ATP ~ A*B*time + (1|heart), data = ATP)
  cap("atp_K_ABt", km(fm_atp, "A:B:time"))
  cap("atp_pm_At", pmo(fm_atp, "A:time", atvar = "time", adj = "BH"))

  data("Drug", package = "predictmeans", envir = environment())
  fm_drug <- glm(cbind(R, N-R) ~ Drug + log10Dose, data = Drug, family = binomial)
  cap("drug_K", km(fm_drug, "Drug"))
  cap("drug_pm", pmo(fm_drug, "Drug", pairwise = TRUE))
  cap("drug_cov", cvm(fm_drug, "Drug", covariate = "log10Dose", trans = binomial()$linkinv)$pltdf)

  data("Clinical", package = "predictmeans", envir = environment())
  fm_clin <- glmer(cbind(Favorable, Total-Favorable) ~ Treatment + (1|Clinic), family = binomial, data = Clinical)
  cap("clin_K", km(fm_clin, "Treatment"))
  cap("clin_pm", pmo(fm_clin, "Treatment", pairwise = TRUE))
  fm_tmb <- glmmTMB::glmmTMB(cbind(Favorable, Total-Favorable) ~ Treatment + (1|Clinic), family = binomial, data = Clinical)
  cap("clin_tmb_K", km(fm_tmb, "Treatment"))

  set.seed(2)
  d <- expand.grid(A = factor(c("a1","a2","a3")), B = factor(c("b1","b2")), rep = 1:6); n <- nrow(d)
  d$x <- runif(n, 1, 10); d$z <- rnorm(n); d$y <- as.numeric(d$A) + d$x + d$z + rnorm(n)
  fm_ax <- lm(y ~ A*x + B + z, d)
  cap("sim_K_cov50", km(fm_ax, "A", "x"))
  cap("sim_K_covnum", km(fm_ax, "A:B", c(2, 0.5)))
  cap("sim_K_covlist", km(fm_ax, "A", c("x","z"), covariateV = list(c(2,5), c(-1,1))))
  cap("sim_K_x_only", km(fm_ax, "x", "x", covariateV = c(2,4,6)))

  data("WarsawApts", package = "HRW", envir = environment())
  W <- WarsawApts; W$district <- factor(W$district)
  sm1 <- suppressWarnings(semireg(areaPerMzloty ~ construction.date,
                 smoothZ = list(grp = smZ(construction.date, k = 25)), data = W, prt = FALSE))
  # 1.1.1 used qnorm(1 - level) for the limits (90% with the default 0.05);
  # 2.0 uses qnorm(1 - level/2), so level = 0.1 gives the limits 1.1.1 gave
  cap("semi_1", semipred(sm1, covariate = "construction.date", level = 0.1, prt = FALSE)$pred_df)
  sm2 <- suppressWarnings(semireg(areaPerMzloty ~ construction.date*district,
                 smoothZ = list(grp = smZ(construction.date, k = 15, by = district, group = TRUE)),
                 data = W, prt = FALSE))
  cap("semi_2_grp", semipred(sm2, "district", "construction.date", level = 0.1, prt = FALSE)$pred_df)
  cap("semi_2_contr", semipred(sm2, "district", "construction.date", contr = c(1,2), level = 0.1, prt = FALSE)$pred_df)
  sm3 <- suppressWarnings(semireg(areaPerMzloty ~ district + surface,
                 smoothZ = list(gc = smZ(construction.date, k = 10)), data = W, prt = FALSE))
  cap("semi_3_smonly", semipred(sm3, covariate = "construction.date", level = 0.1, prt = FALSE)$pred_df)
  cap("semi_3_smonly_grp", semipred(sm3, "district", "construction.date", level = 0.1, prt = FALSE)$pred_df)

  ## Permutation tests, R2, variance components and contrasts
  perm_lmer <- suppressWarnings(permmodels(fm_lmer, method = "raw", nperm = 19, prt = FALSE, ncore = 1, seed = 1))
  cap("perm_lmer", perm_lmer)
  cap("perm_lmer_III", permmodels(fm_lmer, method = "raw", nperm = 9, type = "III", prt = FALSE, ncore = 1, seed = 2)[c("COEFFICENTS", "ANOVA")])
  cap("perm_lme", permmodels(fm_lme, method = "raw", nperm = 9, prt = FALSE, ncore = 1, seed = 3)[c("COEFFICENTS", "ANOVA")])
  cap("perm_lm", permmodels(fm_lm, method = "raw", nperm = 19, prt = FALSE, ncore = 1, seed = 4)[c("COEFFICENTS", "ANOVA")])
  cap("perm_lm_II", permmodels(fm_lm, method = "raw", nperm = 9, type = "II", prt = FALSE, ncore = 1, seed = 5)[c("COEFFICENTS", "ANOVA")])
  cap("perm_pm_lmer", pmo(fm_lmer, "nitro", permlist = perm_lmer, adj = "BH", ncore = 1))
  fm_lmer0 <- lmerTest::lmer(yield ~ nitro + Variety + (1|Block/Whole_plot), data = Oats)
  cap("permlmer", permlmer(fm_lmer0, fm_lmer, nperm = 19, ncore = 1, seed = 6))
  cap("r2_lmer", R2_glmm(fm_lmer))
  cap("r2_lmer_od", R2_glmm(fm_lmer, over_disp = TRUE))
  cap("r2_glmer", R2_glmm(fm_clin))
  cap("r2_tmb", R2_glmm(fm_tmb))
  cap("varcomp_lmer", varcomp(fm_lmer))
  cap("varcomp_lme", varcomp(fm_lme))
  cm <- rbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
  cap("contr_lmer", do.call(contrastmeans, c(list(fm_lmer, "nitro", ctrmatrix = cm), if ("vcov_adjust" %in% names(formals(contrastmeans))) list(vcov_adjust = FALSE))))
  cap("contr_lme", do.call(contrastmeans, c(list(fm_lme, "nitro", ctrmatrix = cm, adj = "BH"), old_df(contrastmeans, fm_lme))))
  R
}

# Drop attributes and row names that are not part of the frozen 1.1.1 output.
clean_reference <- function(x) {
  # R-devel's xtabs() tables carry a "call" attribute; it is not a result
  if (!is.null(attr(x, "call")) && !is.list(x)) attr(x, "call") <- NULL
  if (is.list(x) && !is.data.frame(x)) return(lapply(x, clean_reference))
  # percentages printed as text (R2_glmm) are compared as numbers
  if (is.character(x) && length(x) && all(grepl("^-?[0-9.]+%$", x))) return(as.numeric(sub("%", "", x)))
  if (is.matrix(x)) attr(x, "estimable") <- NULL
  if (is.data.frame(x)) {
    rownames(x) <- NULL
    # lme4 >= 2.0 names the anova column "-2*log(L)" instead of "deviance"
    names(x)[names(x) == "-2*log(L)"] <- "deviance"
    for (v in names(x)) if (is.matrix(x[[v]])) x[[v]] <- as.vector(x[[v]])
  }
  x
}
