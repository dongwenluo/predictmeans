## Calibration of bootmodels() (parametric bootstrap under the global null
## model): type I error of its p-values for a null term when other terms
## have effects (partial null), in LMMs and GLMMs, and under overdispersion.
## Methods: boot (bootmodels, type II Wald chi-square), wald (the same
## statistic against chi-square), kr (Kenward-Roger F; LMM), lrt
## (likelihood-ratio test; GLMM), boot_pair (pair A1 - A2 from the bootstrap
## permlist; split plot). Uses dev/sim-perm-lib.R. Results: dev/sim-boot/.
##
## Run: Rscript dev/sim-boot.R [R] [nsim] [cells|all] [cores]

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500
nperm <- if (length(args) >= 2) as.integer(args[2]) else 199       # bootstrap samples
cores <- if (length(args) >= 4) as.integer(args[4]) else 1
source("dev/sim-perm-lib.R")

description <- list(
  bt_sp = "split plot, 6 blocks x 4 whole plots (A) x 3 (B); A null, B (0,1,2); y ~ A*B + (1|block) + (1|block:A); test A",
  bt_sp_global = "as bt_sp with B null too (the global null bootmodels simulates from)",
  bt_pois = "Poisson GLMM y ~ x + A + (1|block), 8 blocks x 3 A x 2, x lognormal larger for A3, slope 0.15; A null",
  bt_binom3 = "binomial (10 trials) GLMM, 3 blocks, x as bt_pois, logit -0.5 + 0.3 x; A null",
  bt_nb = "negative binomial (size 2) counts fitted as Poisson GLMM, as bt_pois; A null")

boot_p <- function(m, term) {
  bt <- bootmodels(m, nsim = nperm, type = "II", prt = FALSE, ncore = 1)
  list(p = as.numeric(bt$ANOVA[term, "Boot_p_value"]), wald = as.numeric(bt$ANOVA[term, "Pr(>Chisq)"]), bt = bt)
}

one_set <- function(cell) {
  if (cell %in% c("bt_sp", "bt_sp_global")) {
    d <- sp_design()
    d <- d[order(d$block, d$wp, d$sub), ]
    bb <- if (cell == "bt_sp") c(0, 1, 2) else c(0, 0, 0)
    d$y <- bb[d$B] + rnorm(6)[d$block] + rnorm(nlevels(d$wpid))[d$wpid] + rnorm(nrow(d))
    m <- lmer(y ~ A * B + (1 | block) + (1 | block:A), data = d, control = ctl)
    b <- boot_p(m, "A")
    kr <- anova(lmerTest::as_lmerModLmerTest(m), type = "II", ddf = "Kenward-Roger")["A", "Pr(>F)"]
    pr <- predictmeans(m, "A", permlist = b$bt, plot = FALSE, prtnum = FALSE)$pairs
    p12 <- pr$p.value[pr$level1 == "1" & pr$level2 == "2"]
    return(c(boot = b$p, wald = b$wald, kr = kr, boot_pair = p12))
  }
  if (cell == "bt_binom3") {
    d <- lev_design(3); d <- d[rep(seq_len(nrow(d)), each = 2), ]
    d <- d[order(d$block, d$A), ]
    xc <- d$x - mean(d$x)
    d$r <- rbinom(nrow(d), 10, plogis(-0.5 + 0.3 * xc + 0.5 * rnorm(3)[d$block]))
    m <- glmer(cbind(r, 10 - r) ~ x + A + (1 | block), family = binomial, data = d, control = gctl)
    m0 <- glmer(cbind(r, 10 - r) ~ x + (1 | block), family = binomial, data = d, control = gctl)
  } else {
    d <- lev_design(8)
    d <- d[order(d$block, d$A, d$rep), ]
    xc <- d$x - mean(d$x)
    mu <- exp(1 + 0.15 * xc + 0.5 * rnorm(8)[d$block])
    d$y <- if (cell == "bt_nb") rnbinom(nrow(d), size = 2, mu = mu) else rpois(nrow(d), mu)
    m <- glmer(y ~ x + A + (1 | block), family = poisson, data = d, control = gctl)
    m0 <- glmer(y ~ x + (1 | block), family = poisson, data = d, control = gctl)
  }
  b <- boot_p(m, "A")
  lrt <- anova(m0, m)[2, "Pr(>Chisq)"]
  c(boot = b$p, wald = b$wald, lrt = lrt)
}

cells <- if (length(args) >= 3 && args[3] != "all") strsplit(args[3], ",")[[1]] else names(description)
run_cells(cells, one_set, R, nperm, cores, "dev/sim-boot", description)
cat("SIM DONE\n")
