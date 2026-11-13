## Which denominator df are reliable? Type I error at the 5% level of t tests
## of single contrasts and F tests of terms, all true nulls, for every
## df_method of predictmeans, in split plots (lmer, lme) and repeated
## measures (gls, lme), balanced and with missing values. A two-sided t test
## rejecting 5% is a 95% confidence interval covering 95%.
##
## Methods (columns <contrast>.<method>):
##   KR     Kenward-Roger (expected information, pbkrtest), adjusted SEs
##   KRSAS  Kenward-Roger with the observed information (SAS DDFM=KR), adjusted SEs
##   Satt   Satterthwaite, model SEs
##   cont   containment (lmer; SAS DDFM=CONTAIN), model SEs
##   BW     between-within (SAS DDFM=BETWITHIN), model SEs
##   nlme   lme: nlme's df of the term (predictmeans' current lme default)
##   nlmeC  lme: the smallest nlme df of the coefficients each contrast uses
##   resid  gls: residual df (predictmeans' current gls default)
##   z      normal (infinite df)
## Contrasts: see 'description'. Uses dev/sim-perm-lib.R. Results: dev/sim-df/.
##
## Run: Rscript dev/sim-df.R [R] [cells|all] [cores]

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
cores <- if (length(args) >= 3) as.integer(args[3]) else 1
source("dev/sim-perm-lib.R")
inference <- utils::getFromNamespace("pm_inference", "predictmeans")

description <- list(
  sp_bal = "split plot 6 blocks x 4 A (whole plots) x 3 B; var block 1, whole plot 1, residual 1; lmer and lme",
  sp_wp_small = "as sp_bal, whole-plot variance 0.1",
  sp_miss_wp = "as sp_bal, the whole plot (block 1, A1) missing",
  sp_miss_15 = "as sp_bal, 15% of the observations missing at random",
  sp_4blk_miss = "split plot 4 blocks x 4 A x 3 B, 15% missing at random",
  rm_ar1 = "repeated measures, 2 groups x 6 subjects x 4 times, AR(1) rho 0.6, var 1; gls corAR1",
  rm_un = "as rm_ar1 with SDs 1, 1.3, 1.6, 2 and correlation 0.7^lag; gls unstructured (corSymm + varIdent)",
  rm_un_miss = "as rm_un, 10% of the observations missing at random",
  rm_un_12 = "as rm_un with 12 subjects per group",
  rm_cs_lme = "2 groups x 6 subjects x 4 times, subject variance 1, residual 1; lme random intercept, 10% missing")

## contrasts over the means of 'term' as coefficient rows
cells_K <- function(m, term) Kmatrix(m, term)$K

## p-values of the t tests (rows of L) and F tests (list of L) for one method
p_method <- function(m, method, Lt, LF, terms_t, terms_F) {
  dm <- switch(method, KR = "Kenward-Roger", KRSAS = "KR-SAS", Satt = "Satterthwaite", cont = "containment",
               BW = "between-within", "auto")
  inf <- inference(m, dm)
  b <- inf$mp$coef
  V <- inf$mp$vcov
  est <- drop(Lt %*% b)
  se <- sqrt(rowSums((Lt %*% V) * Lt))
  df <- vapply(seq_len(nrow(Lt)), function(i) inf$df(Lt[i, , drop = FALSE], terms_t[i]), numeric(1))
  if (method == "z") df[] <- Inf
  pt_ <- 2 * stats::pt(-abs(est/se), df)
  pF <- vapply(seq_along(LF), function(j) {
    f <- inf$ftest(LF[[j]], terms_F[j])
    if (method == "z") stats::pchisq(f[3] * f[1], f[1], lower.tail = FALSE) else f[4]
  }, numeric(1))
  c(pt_, pF)
}

## nlme's df of the term, and the smallest df of the coefficients each row uses
p_nlme <- function(me, Lt, LF, terms_t, terms_F) {
  b <- nlme::fixef(me)
  V <- as.matrix(stats::vcov(me))
  tdf <- stats::terms(me$fixDF)
  cdf <- me$fixDF$X
  est <- drop(Lt %*% b)
  se <- sqrt(rowSums((Lt %*% V) * Lt))
  dfT <- unname(tdf[terms_t])
  dfC <- apply(Lt, 1, function(l) min(cdf[abs(l) > 1e-10]))
  Ft <- function(L, ddf) {
    q <- qr(t(L))$rank
    L <- L[qr(t(L))$pivot[seq_len(q)], , drop = FALSE]
    bb <- L %*% b
    Fv <- as.numeric(t(bb) %*% solve(L %*% V %*% t(L), bb))/q
    stats::pf(Fv, q, ddf, lower.tail = FALSE)
  }
  pT <- c(2 * stats::pt(-abs(est/se), dfT), vapply(seq_along(LF), function(j) Ft(LF[[j]], tdf[terms_F[j]]), numeric(1)))
  pC <- c(2 * stats::pt(-abs(est/se), dfC),
          vapply(seq_along(LF), function(j) Ft(LF[[j]], min(cdf[colSums(abs(LF[[j]])) > 1e-10])), numeric(1)))
  list(nlme = pT, nlmeC = pC)
}

name_p <- function(p, labs, method) stats::setNames(p, paste(labs, method, sep = "."))

one_set <- function(cell) {
  if (startsWith(cell, "sp_")) {
    nblock <- if (cell == "sp_4blk_miss") 4 else 6
    d <- sp_design(nsub = 3, nblock = nblock)
    s_wp <- if (cell == "sp_wp_small") sqrt(0.1) else 1
    d$y <- rnorm(nblock)[d$block] + s_wp * rnorm(nlevels(d$wpid))[d$wpid] + rnorm(nrow(d))
    if (cell == "sp_miss_wp") d <- d[!(d$block == "1" & d$A == "1"), ]
    if (cell %in% c("sp_miss_15", "sp_4blk_miss")) d <- d[sort(sample(nrow(d), round(0.85 * nrow(d)))), ]
    d <- droplevels(d)
    m <- lmer(y ~ A * B + (1 | block) + (1 | block:A), data = d, control = ctl)
    me <- nlme::lme(y ~ A * B, random = ~ 1 | block/A, data = d)
    KA <- cells_K(m, "A"); KB <- cells_K(m, "B"); KAB <- cells_K(m, "A:B")
    # A1-A2 (whole plot), B1-B2 (subplot), A1B1-A2B1 (whole plot within B1), A1B1-A1B2 (subplot within A1)
    Lt <- rbind(KA[1, ] - KA[2, ], KB[1, ] - KB[2, ], KAB["1:1", ] - KAB["2:1", ], KAB["1:1", ] - KAB["1:2", ])
    LF <- list(A = contrast_family(Kmatrix(m, "A"), "consec") %*% KA,
               AB = t(kronecker(stats::contr.helmert(4), stats::contr.helmert(3))) %*% KAB)
    tt <- c("A", "B", "A:B", "A:B")
    tF <- c("A", "A:B")
    labs <- c("tA", "tB", "tAinB", "tBinA", "FA", "FAB")
    meths <- c("KR", "KRSAS", "Satt", "cont", "BW", "z")
    out <- unlist(lapply(meths, function(k) name_p(p_method(m, k, Lt, LF, tt, tF), labs, k)))
    pn <- p_nlme(me, Lt, LF, tt, tF)
    return(c(out, name_p(pn$nlme, labs, "nlme"), name_p(pn$nlmeC, labs, "nlmeC")))
  }
  ## repeated measures
  nsub <- if (cell == "rm_un_12") 12 else 6
  d <- expand.grid(Time = factor(1:4), Subj = factor(seq_len(2 * nsub)))
  d$Group <- factor(ifelse(as.integer(d$Subj) <= nsub, "A", "B"))
  if (cell == "rm_cs_lme") {
    d$y <- rnorm(2 * nsub)[d$Subj] + rnorm(nrow(d))
  } else {
    sds <- if (cell == "rm_ar1") rep(1, 4) else c(1, 1.3, 1.6, 2)
    rho <- if (cell == "rm_ar1") 0.6 else 0.7
    S <- outer(1:4, 1:4, function(i, j) rho^abs(i - j)) * outer(sds, sds)
    d$y <- as.vector(t(MASS::mvrnorm(2 * nsub, rep(0, 4), S)))
  }
  if (cell %in% c("rm_un_miss", "rm_cs_lme")) d <- d[sort(sample(nrow(d), round(0.9 * nrow(d)))), ]
  d <- droplevels(d)
  d$t <- as.integer(d$Time)
  m <- switch(cell,
    rm_ar1 = nlme::gls(y ~ Group * Time, data = d, correlation = nlme::corAR1(form = ~ t | Subj)),
    rm_cs_lme = nlme::lme(y ~ Group * Time, random = ~ 1 | Subj, data = d),
    nlme::gls(y ~ Group * Time, data = d, correlation = nlme::corSymm(form = ~ t | Subj),
              weights = nlme::varIdent(form = ~ 1 | Time)))
  KG <- cells_K(m, "Group"); KT <- cells_K(m, "Time"); KGT <- cells_K(m, "Group:Time")
  # GA-GB (between), T1-T4 (within), GA-GB at T1 (between at a time), T1-T4 in A (within a group)
  Lt <- rbind(KG[1, ] - KG[2, ], KT[1, ] - KT[4, ], KGT["A:1", ] - KGT["B:1", ], KGT["A:1", ] - KGT["A:4", ])
  LF <- list(T = contrast_family(Kmatrix(m, "Time"), "consec") %*% KT,
             GT = t(kronecker(stats::contr.helmert(2), stats::contr.helmert(4))) %*% KGT)
  tt <- c("Group", "Time", "Group:Time", "Group:Time")
  tF <- c("Time", "Group:Time")
  labs <- c("tG", "tT", "tGatT", "tTinG", "FT", "FGT")
  meths <- c("KR", "KRSAS", "Satt", "BW", "z")
  out <- unlist(lapply(meths, function(k) name_p(p_method(m, k, Lt, LF, tt, tF), labs, k)))
  if (cell == "rm_cs_lme") {
    pn <- p_nlme(m, Lt, LF, tt, tF)
    return(c(out, name_p(pn$nlme, labs, "nlme"), name_p(pn$nlmeC, labs, "nlmeC")))
  }
  c(out, name_p(p_method(m, "auto", Lt, LF, tt, tF), labs, "resid"))
}

cells <- if (length(args) >= 2 && args[2] != "all") strsplit(args[2], ",")[[1]] else names(description)
run_cells(cells, one_set, R, nperm = NA, cores, "dev/sim-df", description)
