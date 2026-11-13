## ML fits of gls: which df? (dev/review-defaults.md, step 3)
## Repeated measures, 2 groups x 6 subjects, unstructured (3 times) or AR(1)
## (4 times), no effects; gls fitted by ML. Type I error at 5% of t tests
## (Group A - B at the last time: between subjects; Time 2 - 1 in group A:
## within) and F tests (Time, Group x Time), with the residual df (the old
## default) or Satterthwaite from the ML log-likelihood (the new default).
## 1000 data sets per design; a true 5% lies in 3.2-6.8% (99%).
## Rscript dev/sim-df-ml.R   (writes dev/sim-df-ml/results.md)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
set.seed(20261114)
nsim <- 1000
one <- function(design) {
  nt <- if (design == "UN") 3 else 4
  d <- data.frame(Subj = factor(rep(1:12, each = nt)), Group = factor(rep(c("A", "B"), each = 6 * nt)),
                  Time = factor(rep(seq_len(nt), 12)), t = rep(seq_len(nt), 12))
  if (design == "UN") {
    S <- matrix(c(1, .6, .3, .6, 1.5, .7, .3, .7, 2), 3)
  } else {
    S <- 0.6^abs(outer(1:4, 1:4, "-"))
  }
  d$y <- as.vector(t(MASS::mvrnorm(12, rep(0, nt), S)))
  g <- tryCatch(if (design == "UN") {
    nlme::gls(y ~ Group * Time, data = d, method = "ML", correlation = nlme::corSymm(form = ~ t | Subj),
              weights = nlme::varIdent(form = ~ 1 | Time))
  } else nlme::gls(y ~ Group * Time, data = d, method = "ML", correlation = nlme::corAR1(form = ~ t | Subj)),
  error = function(e) NULL)
  if (is.null(g)) return(rep(NA, 8))
  K <- Kmatrix(g, "Group:Time")$K
  lev <- rownames(K)
  tB <- K[paste0("A:", nt), ] - K[paste0("B:", nt), ]
  tW <- K["A:2", ] - K["A:1", ]
  last <- cbind(diag(nt - 1), -1)
  FT <- kronecker(rbind(c(.5, .5)), last) %*% K
  FGT <- kronecker(rbind(c(1, -1)), last) %*% K
  out <- c()
  for (m in c("containment", "auto")) {
    inf <- suppressWarnings(suppressMessages(pm_inference(g, m, vcov_adjust = FALSE)))
    b <- inf$mp$coef; V <- inf$mp$vcov
    L <- rbind(tB, tW)[, names(b)]
    tt <- drop(L %*% b) / sqrt(diag(L %*% V %*% t(L)))
    df <- inf$df(L, "Group:Time")
    pf <- c(inf$ftest(FT[, names(b)], "Time")[4], inf$ftest(FGT[, names(b)], "Group:Time")[4])
    out <- c(out, 2 * pt(-abs(tt), df), pf)
  }
  out
}
res <- list()
for (design in c("UN", "AR1")) {
  p <- t(replicate(nsim, one(design)))
  colnames(p) <- paste(rep(c("resid", "Satt"), each = 4), c("t Group", "t Time", "F Time", "F GxT"))
  r <- round(100 * colMeans(p < 0.05, na.rm = TRUE), 1)
  res[[design]] <- data.frame(design = design, test = c("t Group", "t Time", "F Time", "F GxT"),
                              residual = r[1:4], Satterthwaite = r[5:8], n = sum(!is.na(p[, 1])))
}
out <- do.call(rbind, res)
dir.create("dev/sim-df-ml", showWarnings = FALSE)
writeLines(c("# ML gls fits: % of true nulls rejected at 5% (1000 data sets per design)", "",
             "| covariance | test | residual df (old) | Satterthwaite, ML (new) | fits |", "|---|---|---|---|---|",
             apply(out, 1, function(r) paste0("| ", paste(trimws(r), collapse = " | "), " |"))),
           "dev/sim-df-ml/results.md")
print(out, row.names = FALSE)
