## predictmeans against SAS PROC MIXED for AR(1) repeated measures with 5 of 48
## values missing (dev/sas/pm-reference-ar1.sas, run in SAS OnDemand on
## 2026-11-13; parsed results in _data/sas/proc-mixed-ar1.csv.gz). AR(1) is
## not linear in its parameter: SAS's default DDFM=KR adds second-derivative
## terms to the adjusted covariance, KR(FIRSTORDER) leaves them out (as
## predictmeans, which is then free of the parameterisation), KR2 is Kenward
## and Roger (2009). All three give the same df.

ar1_ref <- function() read.csv(test_path("_data", "sas", "proc-mixed-ar1.csv.gz"), colClasses = "character")
ar1_fit <- function() {
  d <- read.table(test_path("_data", "sas", "ar1.txt"), header = TRUE)
  d$Subj <- factor(d$Subj); d$Group <- factor(d$Group); d$t <- d$Time; d$Time <- factor(d$Time)
  nlme::gls(Y ~ Group * Time, data = d, correlation = nlme::corAR1(form = ~ t | Subj))
}
## within one unit of SAS's last printed digit, or 0.1%
ar1_agree <- function(r, s) {
  dec <- ifelse(grepl("\\.", s), nchar(sub("^[^.]*\\.", "", s)), 0)
  v <- as.numeric(s)
  !is.na(r) & (abs(r - v) <= 10^(-dec) * (1 + 1e-6) + 1e-9 | abs(r - v) <= 1e-3 * abs(v))
}
## our values for one SAS run, next to SAS's
ar1_compare <- function(m, S, df_method, se = TRUE) {
  out <- list()
  add <- function(what, key, s, r) out[[length(out) + 1]] <<- data.frame(what = what, key = key, sas = s, r = r)
  for (eff in c("Group", "Time", "Group*Time")) {
    term <- sub("\\*", ":", eff)
    p <- predictmeans(m, term, df_method = df_method, plot = FALSE, prtnum = FALSE)
    L <- S[S$table == "lsmeans" & S$effect == eff, ]
    key <- if (eff == "Group*Time") paste(L$level1, L$level2, sep = ":") else L$level1
    vars <- strsplit(term, ":")[[1]]
    i <- match(key, do.call(paste, c(lapply(p$means[vars], as.character), sep = ":")))
    add("mean", key, L$est, p$means$estimate[i])
    if (se) add("mean se", key, L$se, p$means$std.error[i])
    add("mean df", key, L$df, p$means$df[i])
    D <- S[S$table == "diffs" & S$effect == eff, ]
    k1 <- if (eff == "Group*Time") paste(D$level1, D$level2, sep = ":") else D$level1
    k2 <- if (eff == "Group*Time") paste(D$level3, D$level4, sep = ":") else D$level2
    j <- match(paste(k1, k2), paste(p$pairs$level1, p$pairs$level2))
    if (se) add("diff se", paste(k1, k2), D$se, p$pairs$std.error[j])
    add("diff df", paste(k1, k2), D$df, p$pairs$df[j])
  }
  # Type 3 F tests (pm_joint() uses SAS's basis, each level against the last)
  last <- rbind(c(1, 0, 0, -1), c(0, 1, 0, -1), c(0, 0, 1, -1))
  jt <- rbind(pm_joint(m, "Group", df_method = df_method), pm_joint(m, "Time", df_method = df_method),
              pm_joint(m, "Group:Time", df_method = df_method, ctrmatrix = kronecker(rbind(c(1, -1)), last)))
  T3 <- S[S$table == "tests3", ]
  i <- match(T3$effect, c("Group", "Time", "Group*Time"))
  if (se) add("F", T3$effect, T3$F, jt$F.ratio[i])
  add("den df", T3$effect, T3$dendf, jt$df2[i])
  do.call(rbind, out)
}
ar1_check <- function(res, run) {
  res$ok <- ar1_agree(res$r, res$sas)
  bad <- res[!res$ok, ]
  expect(nrow(bad) == 0, paste0(run, ": ", nrow(bad), " values differ from SAS:\n",
                                paste(utils::capture.output(print(utils::head(bad, 15), row.names = FALSE)), collapse = "\n")))
}

test_that("Kenward-Roger equals SAS DDFM=KR(FIRSTORDER) for AR(1): means, differences, SEs, df, F tests", {
  skip_on_cran()
  m <- ar1_fit()
  ref <- ar1_ref()
  S <- ref[ref$ddfm == "KR(FIRSTORDER)", ]
  cp <- S[S$table == "covparms", ]
  expect_equal(unname(coef(m$modelStruct$corStruct, unconstrained = FALSE)), as.numeric(cp$est[cp$effect == "AR(1)"]),
               tolerance = 2e-4)
  expect_equal(m$sigma^2, as.numeric(cp$est[cp$effect == "Residual"]), tolerance = 2e-4)
  res <- ar1_compare(m, S, "Kenward-Roger")
  expect_gt(nrow(res), 100)
  ar1_check(res, "KR(FIRSTORDER)")
  # SAS's default KR and KR2 give the same df, other standard errors (second-derivative terms)
  for (dd in c("KR", "KR2")) ar1_check(ar1_compare(m, ref[ref$ddfm == dd, ], "Kenward-Roger", se = FALSE), dd)
  kr <- ref[ref$ddfm == "KR" & ref$table == "lsmeans" & ref$effect == "Group", ]
  ours <- predictmeans(m, "Group", plot = FALSE, prtnum = FALSE)$means$std.error
  expect_true(all(ours > as.numeric(kr$se)))                     # 0.3538 against 0.3474
})

test_that("Satterthwaite and between-within equal SAS for AR(1) with missing values", {
  skip_on_cran()
  m <- ar1_fit()
  ref <- ar1_ref()
  ar1_check(ar1_compare(m, ref[ref$ddfm == "SATTERTHWAITE", ], "Satterthwaite"), "SATTERTHWAITE")
  ar1_check(ar1_compare(m, ref[ref$ddfm == "BETWITHIN", ], "between-within"), "BETWITHIN")
})
