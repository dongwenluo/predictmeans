## Follow-up of dev/sim-df.R: Kenward-Roger with the expected information
## (pbkrtest; df_method = "Kenward-Roger") or the observed information (SAS
## DDFM=KR; "KR-SAS") when values are missing. With unstructured repeated
## measures and 10% missing the expected information gave 8.5% for F tests,
## the observed 6% (dev/sim-df/results.md, rm_un_miss). Kenward and
## Molenberghs (1998) show the expected information is not valid when data
## are missing (at random); does the difference matter for lmer too?
## Methods and contrasts as dev/sim-df.R. Results: dev/sim-df-kr/.
##
## Run: Rscript dev/sim-df-kr.R [R] [cells|all] [cores]

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
cores <- if (length(args) >= 3) as.integer(args[3]) else 1
source("dev/sim-perm-lib.R")
inference <- utils::getFromNamespace("pm_inference", "predictmeans")

description <- list(
  sp_miss_30 = "split plot 6 blocks x 4 A x 3 B (variances 1, 1, 1), 30% of the observations missing; lmer",
  rs_lmer_miss = "2 groups x 8 subjects x 4 times; random intercept (var 1) and slope in time (var 0.2), residual 1; 20% missing; lmer (t | Subj)",
  cs_lmer_miss = "2 groups x 6 subjects x 4 times; subject var 1, residual 1; 20% missing; lmer (1 | Subj)",
  rm_un_miss20 = "2 groups x 12 subjects x 4 times; SDs 1, 1.3, 1.6, 2, correlation 0.7^lag; 20% missing; gls unstructured",
  rm_ar1_miss20 = "2 groups x 6 subjects x 6 times; AR(1) rho 0.6, SDs 1 to 2; 20% missing; gls corAR1 + varIdent")

p_method <- function(m, method, Lt, LF, terms_t, terms_F) {
  dm <- switch(method, KR = "Kenward-Roger", KRSAS = "KR-SAS", Satt = "Satterthwaite", cont = "containment",
               BW = "between-within", "auto")
  inf <- inference(m, dm)
  b <- inf$mp$coef
  V <- inf$mp$vcov
  est <- drop(Lt %*% b)
  se <- sqrt(rowSums((Lt %*% V) * Lt))
  df <- vapply(seq_len(nrow(Lt)), function(i) inf$df(Lt[i, , drop = FALSE], terms_t[i]), numeric(1))
  pt_ <- 2 * stats::pt(-abs(est/se), df)
  pF <- vapply(seq_along(LF), function(j) inf$ftest(LF[[j]], terms_F[j])[4], numeric(1))
  c(pt_, pF)
}
name_p <- function(p, labs, method) stats::setNames(p, paste(labs, method, sep = "."))

one_set <- function(cell) {
  if (cell == "sp_miss_30") {
    d <- sp_design(nsub = 3, nblock = 6)
    d$y <- rnorm(6)[d$block] + rnorm(nlevels(d$wpid))[d$wpid] + rnorm(nrow(d))
    d <- droplevels(d[sort(sample(nrow(d), round(0.7 * nrow(d)))), ])
    m <- lmer(y ~ A * B + (1 | block) + (1 | block:A), data = d, control = ctl)
    KA <- Kmatrix(m, "A")$K; KB <- Kmatrix(m, "B")$K; KAB <- Kmatrix(m, "A:B")$K
    Lt <- rbind(KA[1, ] - KA[2, ], KB[1, ] - KB[2, ], KAB["1:1", ] - KAB["2:1", ], KAB["1:1", ] - KAB["1:2", ])
    LF <- list(contrast_family(Kmatrix(m, "A"), "consec") %*% KA,
               t(kronecker(stats::contr.helmert(4), stats::contr.helmert(3))) %*% KAB)
    tt <- c("A", "B", "A:B", "A:B"); tF <- c("A", "A:B")
    labs <- c("tA", "tB", "tAinB", "tBinA", "FA", "FAB")
    meths <- c("KR", "KRSAS", "Satt", "cont")
  } else {
    nsub <- switch(cell, rs_lmer_miss = 8, rm_un_miss20 = 12, 6)
    nt <- if (cell == "rm_ar1_miss20") 6 else 4
    d <- expand.grid(Time = factor(seq_len(nt)), Subj = factor(seq_len(2 * nsub)))
    d$Group <- factor(ifelse(as.integer(d$Subj) <= nsub, "A", "B"))
    d$t <- as.integer(d$Time)
    ns <- 2 * nsub
    if (cell == "rs_lmer_miss") {
      d$y <- rnorm(ns)[d$Subj] + sqrt(0.2) * rnorm(ns)[d$Subj] * (d$t - 2.5) + rnorm(nrow(d))
    } else if (cell == "cs_lmer_miss") {
      d$y <- rnorm(ns)[d$Subj] + rnorm(nrow(d))
    } else {
      sds <- seq(1, 2, length.out = nt)
      rho <- if (cell == "rm_ar1_miss20") 0.6 else 0.7
      S <- outer(seq_len(nt), seq_len(nt), function(i, j) rho^abs(i - j)) * outer(sds, sds)
      d$y <- as.vector(t(MASS::mvrnorm(ns, rep(0, nt), S)))
    }
    d <- droplevels(d[sort(sample(nrow(d), round(0.8 * nrow(d)))), ])
    m <- switch(cell,
      rs_lmer_miss = lmer(y ~ Group * Time + (t | Subj), data = d, control = ctl),
      cs_lmer_miss = lmer(y ~ Group * Time + (1 | Subj), data = d, control = ctl),
      rm_ar1_miss20 = nlme::gls(y ~ Group * Time, data = d, correlation = nlme::corAR1(form = ~ t | Subj),
                                weights = nlme::varIdent(form = ~ 1 | Time)),
      nlme::gls(y ~ Group * Time, data = d, correlation = nlme::corSymm(form = ~ t | Subj),
                weights = nlme::varIdent(form = ~ 1 | Time)))
    KG <- Kmatrix(m, "Group")$K; KT <- Kmatrix(m, "Time")$K; KGT <- Kmatrix(m, "Group:Time")$K
    last <- as.character(nt)
    Lt <- rbind(KG[1, ] - KG[2, ], KT[1, ] - KT[nt, ], KGT["A:1", ] - KGT["B:1", ],
                KGT["A:1", ] - KGT[paste0("A:", last), ])
    LF <- list(contrast_family(Kmatrix(m, "Time"), "consec") %*% KT,
               t(kronecker(stats::contr.helmert(2), stats::contr.helmert(nt))) %*% KGT)
    tt <- c("Group", "Time", "Group:Time", "Group:Time"); tF <- c("Time", "Group:Time")
    labs <- c("tG", "tT", "tGatT", "tTinG", "FT", "FGT")
    meths <- c("KR", "KRSAS", "Satt", "BW")
  }
  unlist(lapply(meths, function(k) name_p(p_method(m, k, Lt, LF, tt, tF), labs, k)))
}

cells <- if (length(args) >= 2 && args[2] != "all") strsplit(args[2], ",")[[1]] else names(description)
run_cells(cells, one_set, R, nperm = NA, cores, "dev/sim-df-kr", description)
