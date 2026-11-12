## Study 1b, part 2: tests of a fixed term by permutation with the package's
## own functions (plan §10 step 1b; chooses the defaults of step 9).
## Uses dev/sim-perm-lib.R. Results: dev/sim-perm-fixed/.
##
## Methods (p-value of the tested term A in each data set):
##   raw       permmodels(), the 1.x scheme (lm/glm: response permuted;
##             lmer/lme: fitted values and decorrelated residuals together;
##             refits)
##   fl        permmodels(method = "freedman-lane", refit = FALSE): V fixed at
##             the fit, decorrelated least squares (GLM/GLMM: on the working
##             scale), no refits
##   flr       permmodels(method = "freedman-lane") for lmer/lme: the model
##             without the term fitted, its decorrelated residuals permuted,
##             the full model refitted to each permuted data set; for GLMMs
##             the same with the working linear mixed model
##   (pairs cell) flr: pairs by Freedman-Lane with refits from the model with
##             all the factor's means equal (predictmeans(permlist = )), kr:
##             predictmeans() t tests
##   permlmer  permlmer(H0, H1) for GLMMs (working scale with weighted lmer
##             refits; quantile residuals and GLMM refits for 0/1 data)
##   kr / lrt / param   the usual test: Kenward-Roger F (lmer), likelihood
##             ratio chi-square (glm, glmer, glmmTMB), nlme's F (lme)
##
## Run: Rscript dev/sim-perm-fixed.R [R] [nperm] [cells|all] [cores]

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
nperm <- if (length(args) >= 2) as.integer(args[2]) else 199
cores <- if (length(args) >= 4) as.integer(args[4]) else 1
source("dev/sim-perm-lib.R")
suppressPackageStartupMessages({ library(nlme); library(glmmTMB) })

pA <- function(tab, term = "A") as.numeric(tab[term, "Perm_p_value"])

lmm_term <- function(d, f, type = "II") {
  m <- lmer(f, data = d, control = ctl)
  raw <- permmodels(m, nperm = nperm, type = type, prt = FALSE, ncore = 1)$ANOVA
  fl <- permmodels(m, nperm = nperm, type = type, prt = FALSE, method = "freedman-lane", refit = FALSE)$ANOVA
  # (refits for the tested term only, to save time)
  flr <- predictmeans:::pm_perm_fl_refit(m, nperm, type, FALSE, NULL, FALSE, NULL, 1, only = "A")$ANOVA
  kr <- anova(lmerTest::as_lmerModLmerTest(m), type = type, ddf = "Kenward-Roger")["A", "Pr(>F)"]
  c(raw = pA(raw), fl = pA(fl), flr = pA(flr), kr = kr)
}

glm_term <- function(d, f, family) {
  m <- glm(f, family = family, data = d)                      # A is the last term
  raw <- permmodels(m, nperm = nperm, type = "I", prt = FALSE, ncore = 1)$ANOVA
  fl <- permmodels(m, nperm = nperm, type = "II", prt = FALSE, method = "freedman-lane")$ANOVA
  lrt <- anova(m, test = "Chisq")["A", "Pr(>Chi)"]
  c(raw = pA(raw), fl = pA(fl), lrt = lrt)
}

glmm_term <- function(d, f0, f1, family, tmb = FALSE) {
  fit <- function(f) if (tmb) glmmTMB(f, family = family, data = d) else glmer(f, family = family, data = d, control = gctl)
  g0 <- fit(f0)
  g1 <- fit(f1)
  fl <- permmodels(g1, nperm = nperm, type = "II", prt = FALSE, method = "freedman-lane", refit = FALSE)$ANOVA
  # refits of the working LMM, for the tested term only (to save time)
  flr <- predictmeans:::pm_perm_fl_refit(g1, nperm, "II", FALSE, NULL, FALSE, NULL, 1, only = "A")$ANOVA
  pl <- permlmer(g0, g1, nperm = nperm, ncore = 1)$`Perm-p`[2]
  df <- length(pm_fixef(g1)) - length(pm_fixef(g0))
  lrt <- stats::pchisq(2 * (as.numeric(logLik(g1)) - as.numeric(logLik(g0))), df, lower.tail = FALSE)
  c(fl = pA(fl), flr = pA(flr), permlmer = pl, lrt = lrt)
}
pm_fixef <- predictmeans:::pm_fixef

## rows: differences of the A means (averaged over B), pairs i < j
pair_rows <- function(m, d) {
  grid <- expand.grid(A = factor(levels(d$A), levels = levels(d$A)), B = factor(levels(d$B), levels = levels(d$B)))
  X <- model.matrix(~ A * B, grid)
  M <- rowsum(X, grid$A) / nlevels(d$B)
  pr <- utils::combn(nlevels(d$A), 2)
  L <- M[pr[2, ], ] - M[pr[1, ], ]
  rownames(L) <- paste0(pr[1, ], "v", pr[2, ])
  L
}

rcbd <- function(nblock = 8, reps = 1) {
  d <- expand.grid(rep = seq_len(reps), A = factor(1:4), block = factor(seq_len(nblock)))
  d[order(d$block, d$A), ]
}

description <- list(
  fx2_lmm_lev = "lmer y ~ x + A + (1|block), 8 blocks x 3 A x 2, x lognormal larger for A3, slope 1; A null",
  fx2_lmm_sp = "split plot, 6 blocks, A (whole plot) null, B (0,1,2); y ~ A*B + (1|block) + (1|block:A); A",
  fx2_lmm_unbal = "lmer y ~ x + A + (1|block), 8 blocks x 4 A x 2 with 20% missing, x normal related to A; A null",
  fx2_lme_varident = "lme y ~ A, random 1|block, varIdent by A (sd 1,1,2,3), 8 blocks x 4 x 2; A null",
  fx2_glm_pois_lev = "glm Poisson y ~ x + A, x lognormal larger for A3, log mean 1 + 0.15 x; A null",
  fx2_glm_binom_lev = "glm binomial (10 trials) y ~ x + A, x as above, logit -0.5 + 0.3 x; A null",
  fx2_glm_nb = "counts negative binomial (size 2) fitted as Poisson glm, x as above; A null",
  fx2_glmm_pois_lev = "Poisson GLMM y ~ x + A + (1|block), 8 blocks, x as above, block sd 0.5; A null",
  fx2_glmm_nb = "negative binomial (size 2) fitted as Poisson GLMM, as fx2_glmm_pois_lev; A null",
  fx2_glmm_binom3 = "binomial (10 trials) GLMM, 3 blocks, x as above; A null",
  fx2_glmm_rare = "binary GLMM, about 5-10% events, A tested with B associated, 8 blocks x 3; A null",
  fx2_tmb_nb2 = "glmmTMB nbinom2 (size 2) y ~ x + A + (1|block), as fx2_glmm_nb but the right family; A null",
  fx2_lmm_lev_power = "as fx2_lmm_lev with A effects (0, 0.5, 1)",
  fx2_glmm_pois_power = "as fx2_glmm_pois_lev with A effects (0, 0.2, 0.4) on the log scale",
  fx2_sp_pairs = "split plot, 6 blocks, A (whole plot) means (0,0,0,2), B (0,1,2); pairs of A means by Freedman-Lane with refits (as pr_sp_A_partial)")

one_set <- function(cell) {
  switch(cell,
    fx2_lmm_lev = , fx2_lmm_lev_power = {
      d <- lev_design(8)
      d <- d[order(d$block, d$A, d$rep), ]
      eff <- if (cell == "fx2_lmm_lev_power") c(0, 0.5, 1)[d$A] else 0
      d$y <- d$x + eff + rnorm(8)[d$block] + rnorm(nrow(d))
      lmm_term(d, y ~ x + A + (1 | block))
    },
    fx2_lmm_sp = {
      d <- sp_design()
      d <- d[order(d$block, d$wp, d$sub), ]
      d$y <- c(0, 1, 2)[d$B] + rnorm(6)[d$block] + rnorm(nlevels(d$wpid))[d$wpid] + rnorm(nrow(d))
      lmm_term(d, y ~ A * B + (1 | block) + (1 | block:A))
    },
    fx2_lmm_unbal = {
      d <- rcbd(8, 2)
      d$x <- rnorm(nrow(d)) + 0.5 * (d$A == "4")
      d$y <- d$x + rnorm(8)[d$block] + rnorm(nrow(d))
      d <- d[sort(sample(nrow(d), round(0.8 * nrow(d)))), ]
      lmm_term(d, y ~ x + A + (1 | block))
    },
    fx2_lme_varident = {
      d <- rcbd(8, 2)
      d$y <- rnorm(8)[d$block] + c(1, 1, 2, 3)[d$A] * rnorm(nrow(d))
      m <- lme(y ~ A, random = ~ 1 | block, weights = varIdent(form = ~ 1 | A), data = d)
      raw <- permmodels(m, nperm = nperm, type = "II", prt = FALSE, ncore = 1)$ANOVA
      fl <- permmodels(m, nperm = nperm, type = "II", prt = FALSE, method = "freedman-lane", refit = FALSE)$ANOVA
      flr <- predictmeans:::pm_perm_fl_refit(m, nperm, "II", FALSE, NULL, FALSE, NULL, 1, only = "A")$ANOVA
      c(raw = pA(raw), fl = pA(fl), flr = pA(flr), param = anova(m)["A", "p-value"])
    },
    fx2_glm_pois_lev = , fx2_glm_binom_lev = , fx2_glm_nb = {
      d <- lev_design(8)
      xc <- d$x - mean(d$x)
      if (cell == "fx2_glm_binom_lev") {
        d$r <- rbinom(nrow(d), 10, plogis(-0.5 + 0.3 * xc))
        glm_term(d, cbind(r, 10 - r) ~ x + A, binomial)
      } else {
        mu <- exp(1 + 0.15 * xc)
        d$y <- if (cell == "fx2_glm_nb") rnbinom(nrow(d), size = 2, mu = mu) else rpois(nrow(d), mu)
        glm_term(d, y ~ x + A, poisson)
      }
    },
    fx2_glmm_pois_lev = , fx2_glmm_nb = , fx2_tmb_nb2 = , fx2_glmm_pois_power = {
      d <- lev_design(8)
      d <- d[order(d$block, d$A, d$rep), ]
      xc <- d$x - mean(d$x)
      eff <- if (cell == "fx2_glmm_pois_power") c(0, 0.2, 0.4)[d$A] else 0
      mu <- exp(1 + 0.15 * xc + eff + 0.5 * rnorm(8)[d$block])
      d$y <- if (cell %in% c("fx2_glmm_nb", "fx2_tmb_nb2")) rnbinom(nrow(d), size = 2, mu = mu) else rpois(nrow(d), mu)
      if (cell == "fx2_tmb_nb2") {
        glmm_term(d, y ~ x + (1 | block), y ~ x + A + (1 | block), nbinom2, tmb = TRUE)
      } else glmm_term(d, y ~ x + (1 | block), y ~ x + A + (1 | block), poisson)
    },
    fx2_sp_pairs = {
      d <- sp_design()
      d <- d[order(d$block, d$wp, d$sub), ]
      d$y <- c(0, 0, 0, 2)[d$A] + c(0, 1, 2)[d$B] + rnorm(6)[d$block] + rnorm(nlevels(d$wpid))[d$wpid] + rnorm(nrow(d))
      m <- lmer(y ~ A * B + (1 | block) + (1 | block:A), data = d, control = ctl)
      perm <- structure(list(scheme = "freedman-lane-refit", nperm = nperm, exact = FALSE,
                             start = attr(predictmeans:::pm_perm_stream(nrow(d)), "start"), idx = NULL,
                             coef_names = names(fixef(m))), class = "pm_perm")
      pr <- as.data.frame(predictmeans(m, "A", plot = FALSE, prtnum = FALSE)$pairs)
      L <- pair_rows(m, d)
      p_flr <- predictmeans:::pm_perm_fl_refit_t(perm, L, m, 1)
      key <- paste0(pr$level1, "v", pr$level2)
      p_kr <- pr$p.value[match(rownames(L), key)]
      nul <- rownames(L) %in% c("1v2", "1v3", "2v3")
      c(setNames(p_flr, paste0("flr.", rownames(L))), setNames(p_kr, paste0("kr.", rownames(L))),
        flr.holm = min(p.adjust(p_flr, "holm")[nul]), kr.holm = min(p.adjust(p_kr, "holm")[nul]))
    },
    fx2_glmm_binom3 = {
      d <- lev_design(3); d <- d[rep(seq_len(nrow(d)), each = 2), ]
      d <- d[order(d$block, d$A), ]
      xc <- d$x - mean(d$x)
      d$r <- rbinom(nrow(d), 10, plogis(-0.5 + 0.3 * xc + 0.5 * rnorm(3)[d$block]))
      glmm_term(d, cbind(r, 10 - r) ~ x + (1 | block), cbind(r, 10 - r) ~ x + A + (1 | block), binomial)
    },
    fx2_glmm_rare = {
      d <- assoc_design(8, reps = 3)
      d$y <- rbinom(nrow(d), 1, plogis(c(-3.5, -2.8, -2.1)[d$B] + 0.5 * rnorm(8)[d$block]))
      glmm_term(d, y ~ B + (1 | block), y ~ B + A + (1 | block), binomial)
    },
    stop("unknown cell ", cell))
}

cells <- if (length(args) >= 3 && args[3] != "all") strsplit(args[3], ",")[[1]] else names(description)
run_cells(cells, one_set, R, nperm, cores, "dev/sim-perm-fixed", description)
cat("SIM DONE\n")
