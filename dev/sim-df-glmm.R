## Denominator df for GLMMs: type I error at the 5% level of Wald t tests of
## single contrasts and Wald F tests of terms (all true nulls) in glmer fits,
## with the df predictmeans uses now (residual, n - parameters) and the
## alternatives: containment (SAS GLIMMIX's default with random effects),
## between-within, nlme's df (the level at which each term varies: 'nlme' by
## term, 'nlmeC' the smallest of the coefficients a contrast uses; computed by
## nlme from the design, which is all they depend on), and none (normal /
## chi-square, as emmeans). GLMMs have no Kenward-Roger or Satterthwaite df here. Uses dev/sim-perm-lib.R.
## Results: dev/sim-df-glmm/.
##
## Run: Rscript dev/sim-df-glmm.R [R] [cells|all] [cores]

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
cores <- if (length(args) >= 3) as.integer(args[3]) else 1
source("dev/sim-perm-lib.R")
cont_df <- utils::getFromNamespace("pm_containment_df", "predictmeans")
bw_df <- utils::getFromNamespace("pm_betwithin_df", "predictmeans")

description <- list(
  gl_pois_sp = "Poisson split plot, 6 blocks x 4 A (whole plots) x 3 B; log mean 2; SD block 0.3, whole plot 0.3",
  gl_pois_sp_low = "as gl_pois_sp with log mean 0.5 (mean count 1.6)",
  gl_pois_sp_4blk = "as gl_pois_sp with 4 blocks",
  gl_binom_clus = "binomial (10 trials), 2 treatments x 4 clusters x 6 units, within-cluster factor W (2 levels); logit -0.5; cluster SD 0.5",
  gl_binom_clus_8 = "as gl_binom_clus with 8 clusters per treatment",
  gl_bern_rcbd = "binary, 8 blocks x 3 treatments x 4 units; logit 0; block SD 0.5",
  gl_pois_crossed = "Poisson, 4 rows x 12 columns crossed (random), treatment T (3 levels) on columns, W (2 levels) within each row x column; log mean 1.5; SD row 0.3, column 0.3; df methods of the package (inner-outer extended to crossed factors)",
  gl_binom_crossed = "binomial (10 trials), 6 rows x 8 columns crossed, T (2 levels) on columns, W within cells; logit -0.5; SD row 0.5, column 0.5")

## crossed grouping factors: the package's own df methods (pm_inference())
inference <- utils::getFromNamespace("pm_inference", "predictmeans")
one_set_crossed <- function(cell) {
  pois <- cell == "gl_pois_crossed"
  nr <- if (pois) 4 else 6; nc <- if (pois) 12 else 8; nT <- if (pois) 3 else 2
  d <- expand.grid(W = factor(1:2), Row = factor(seq_len(nr)), Col = factor(seq_len(nc)))
  d$T <- factor((as.integer(d$Col) - 1) %% nT + 1)
  if (pois) {
    d$y <- rpois(nrow(d), exp(1.5 + 0.3 * rnorm(nr)[d$Row] + 0.3 * rnorm(nc)[d$Col]))
    m <- glmer(y ~ T * W + (1 | Row) + (1 | Col), family = poisson, data = d, control = gctl)
  } else {
    d$r <- rbinom(nrow(d), 10, plogis(-0.5 + 0.5 * rnorm(nr)[d$Row] + 0.5 * rnorm(nc)[d$Col]))
    m <- glmer(cbind(r, 10 - r) ~ T * W + (1 | Row) + (1 | Col), family = binomial, data = d, control = gctl)
  }
  KT <- Kmatrix(m, "T")$K; KW <- Kmatrix(m, "W")$K; KTW <- Kmatrix(m, "T:W")$K
  Lt <- rbind(KT[1, ] - KT[2, ], KW[1, ] - KW[2, ], KTW["1:1", ] - KTW["2:1", ], KTW["1:1", ] - KTW["1:2", ])
  LF <- list(contrast_family(Kmatrix(m, "T"), "consec") %*% KT)
  tt <- c("T", "W", "T:W", "T:W"); labs <- c("tT", "tW", "tTinW", "tWinT", "FT")
  b <- lme4::fixef(m); V <- as.matrix(stats::vcov(m))
  est <- drop(Lt %*% b); se <- sqrt(rowSums((Lt %*% V) * Lt))
  out <- unlist(lapply(c("inner-outer", "residual", "containment", "z"), function(k) {
    inf <- inference(m, if (k == "z") "residual" else k)
    df <- vapply(seq_len(nrow(Lt)), function(i) inf$df(Lt[i, , drop = FALSE], tt[i]), numeric(1))
    f <- inf$ftest(LF[[1]], "T")
    if (k == "z") { df[] <- Inf; f[4] <- stats::pchisq(f[3] * f[1], f[1], lower.tail = FALSE) }
    p <- c(2 * stats::pt(-abs(est/se), df), f[4])
    stats::setNames(p, paste(labs, c("inner-outer" = "io", residual = "resid", containment = "cont", z = "z")[[k]], sep = "."))
  }))
  inf <- inference(m, "inner-outer")
  c(out, df_io_tT = inf$df(Lt[1, , drop = FALSE], "T"), df_io_tW = inf$df(Lt[2, , drop = FALSE], "W"))
}

one_set <- function(cell) {
  if (cell %in% c("gl_pois_crossed", "gl_binom_crossed")) return(one_set_crossed(cell))
  if (startsWith(cell, "gl_pois_sp")) {
    nb <- if (cell == "gl_pois_sp_4blk") 4 else 6
    d <- sp_design(nsub = 3, nblock = nb)
    mu0 <- if (cell == "gl_pois_sp_low") 0.5 else 2
    d$y <- rpois(nrow(d), exp(mu0 + 0.3 * rnorm(nb)[d$block] + 0.3 * rnorm(nlevels(d$wpid))[d$wpid]))
    m <- glmer(y ~ A * B + (1 | block) + (1 | block:A), family = poisson, data = d, control = gctl)
    d$z <- rnorm(nrow(d))
    me <- nlme::lme(z ~ A * B, random = ~ 1 | block/A, data = d, control = nlme::lmeControl(returnObject = TRUE))
    KA <- Kmatrix(m, "A")$K; KB <- Kmatrix(m, "B")$K; KAB <- Kmatrix(m, "A:B")$K
    Lt <- rbind(KA[1, ] - KA[2, ], KB[1, ] - KB[2, ], KAB["1:1", ] - KAB["2:1", ], KAB["1:1", ] - KAB["1:2", ])
    LF <- list(contrast_family(Kmatrix(m, "A"), "consec") %*% KA,
               t(kronecker(stats::contr.helmert(4), stats::contr.helmert(3))) %*% KAB)
    tt <- c("A", "B", "A:B", "A:B"); tF <- c("A", "A:B")
    labs <- c("tA", "tB", "tAinB", "tBinA", "FA", "FAB")
  } else {
    if (cell == "gl_bern_rcbd") {
      d <- expand.grid(rep = 1:4, Trt = factor(1:3), block = factor(1:8))
      d$y <- rbinom(nrow(d), 1, plogis(0.5 * rnorm(8)[d$block]))
      m <- glmer(y ~ Trt + (1 | block), family = binomial, data = d, control = gctl)
      d$z <- rnorm(nrow(d))
      me <- nlme::lme(z ~ Trt, random = ~ 1 | block, data = d, control = nlme::lmeControl(returnObject = TRUE))
      K <- Kmatrix(m, "Trt")$K
      Lt <- rbind(K[1, ] - K[2, ], K[1, ] - K[3, ])
      LF <- list(contrast_family(Kmatrix(m, "Trt"), "consec") %*% K)
      tt <- c("Trt", "Trt"); tF <- "Trt"
      labs <- c("t12", "t13", "FTrt")
    } else {
      nc <- if (cell == "gl_binom_clus_8") 8 else 4
      d <- expand.grid(unit = 1:3, W = factor(1:2), clus = factor(seq_len(2 * nc)))
      d$Trt <- factor(ifelse(as.integer(d$clus) <= nc, "a", "b"))
      d$r <- rbinom(nrow(d), 10, plogis(-0.5 + 0.5 * rnorm(2 * nc)[d$clus]))
      m <- glmer(cbind(r, 10 - r) ~ Trt * W + (1 | clus), family = binomial, data = d, control = gctl)
      d$z <- rnorm(nrow(d))
      me <- nlme::lme(z ~ Trt * W, random = ~ 1 | clus, data = d, control = nlme::lmeControl(returnObject = TRUE))
      KT <- Kmatrix(m, "Trt")$K; KW <- Kmatrix(m, "W")$K; KTW <- Kmatrix(m, "Trt:W")$K
      Lt <- rbind(KT[1, ] - KT[2, ], KW[1, ] - KW[2, ], KTW["a:1", ] - KTW["b:1", ])
      LF <- list(t(kronecker(stats::contr.helmert(2), stats::contr.helmert(2))) %*% KTW)
      tt <- c("Trt", "W", "Trt:W"); tF <- "Trt:W"
      labs <- c("tTrt", "tW", "tTrtinW", "FTW")
    }
  }
  b <- lme4::fixef(m)
  V <- as.matrix(stats::vcov(m))
  est <- drop(Lt %*% b)
  se <- sqrt(rowSums((Lt %*% V) * Lt))
  W <- vapply(LF, function(L) {
    q <- qr(t(L))$rank
    L <- L[qr(t(L))$pivot[seq_len(q)], , drop = FALSE]
    bb <- L %*% b
    c(q, as.numeric(t(bb) %*% solve(L %*% V %*% t(L), bb)))
  }, numeric(2))
  cd <- cont_df(m)
  bd <- bw_df(m)
  tdf <- stats::terms(me$fixDF)
  cdf <- me$fixDF$X
  dfs <- list(resid = rep(summary(m)$AICtab[["df.resid"]], length(tt) + length(tF)),
              cont = unname(c(cd[tt], cd[tF])), BW = unname(c(bd[tt], bd[tF])),
              nlme = unname(c(tdf[tt], tdf[tF])),
              nlmeC = c(apply(Lt, 1, function(l) min(cdf[abs(l) > 1e-10])),
                        vapply(LF, function(L) min(cdf[colSums(abs(L)) > 1e-10]), numeric(1))),
              z = rep(Inf, length(tt) + length(tF)))
  out <- unlist(lapply(names(dfs), function(k) {
    dd <- dfs[[k]]
    nt <- length(tt)
    p <- c(2 * stats::pt(-abs(est/se), dd[seq_len(nt)]),
           stats::pf(W[2, ]/W[1, ], W[1, ], dd[nt + seq_along(tF)], lower.tail = FALSE))
    stats::setNames(p, paste(labs, k, sep = "."))
  }))
  # the df themselves, to report (not p-values: excluded from rejection rates by the summary)
  c(out, stats::setNames(unlist(lapply(dfs, `[`, 1:2)), paste0("df_", rep(names(dfs), each = 2), "_", labs[1:2])))
}

cells <- if (length(args) >= 2 && args[2] != "all") strsplit(args[2], ",")[[1]] else names(description)
run_cells(cells, one_set, R, nperm = NA, cores, "dev/sim-df-glmm", description)
