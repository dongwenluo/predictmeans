## predictmeans against SAS PROC MIXED (reference values: _data/sas/README.md).
## Every variance component, LS-mean, difference and Type 3 test of the 16
## runs must agree with what SAS printed: within one unit of its last digit
## or 0.1% (REML iterations stop at slightly different places; the 1% to 10%
## between "KR-expected" and SAS's KR is still caught).

sas_ref <- function() read.csv(test_path("_data", "sas", "proc-mixed-reference.csv.gz"), colClasses = "character")

sas_oats <- function() {
  d <- as.data.frame(nlme::Oats)
  data.frame(Block = factor(as.character(d$Block)), Variety = factor(gsub(" ", "", d$Variety)),
             Nitro = factor(sprintf("%.1f", d$nitro)), Yield = d$yield)
}
sas_rep <- function() {
  y <- c(20.1, 24.5, 28.2, 18.4, 22.0, 25.9, 22.3, 26.1, 30.4, 19.0, 23.4, 27.1, 21.5, 25.0, 29.0, 17.8, 21.2, 24.8,
         19.5, 21.0, 22.4, 21.0, 22.8, 23.9, 18.2, 19.5, 20.8, 20.4, 21.9, 23.1, 22.1, 23.5, 24.6, 19.0, 20.1, 21.5)
  data.frame(Subj = factor(rep(1:12, each = 3)), Group = factor(rep(c("TrtA", "TrtB"), each = 18)),
             Time = factor(rep(1:3, 12)), Y = y)
}

## agreement with a printed SAS value
sas_agree <- function(r, s) {
  ok <- s != "" & !grepl("<", s)
  dec <- ifelse(grepl("\\.", s), nchar(sub("^[^.]*\\.", "", s)), 0)
  v <- suppressWarnings(as.numeric(s))
  out <- rep(TRUE, length(s))
  out[ok] <- !is.na(r[ok]) & (abs(r[ok] - v[ok]) <= 10^(-dec[ok]) * (1 + 1e-6) + 1e-9 |
                                abs(r[ok] - v[ok]) <= 1e-3 * abs(v[ok]))
  out
}

## the comparisons of one SAS run: data frame (what, key, sas, r)
sas_run <- function(S) {
  dname <- S$data[1]
  ddfm <- S$ddfm[1]
  res <- list()
  add <- function(what, key, s, r) res[[length(res) + 1]] <<- data.frame(what = what, key = key, sas = s, r = r)
  pair_keys <- function(D, two) {
    k1 <- if (two) paste(D$level1, D$level2, sep = ":") else D$level1
    k2 <- if (two) paste(D$level3, D$level4, sep = ":") else D$level2
    list(k1, k2)
  }
  if (dname != "REPEATEDLONG") {
    d <- sas_oats()
    if (dname %in% c("OATS_MISSWP", "OATS_MISS")) d <- d[!(d$Block == "I" & d$Variety == "GoldenRain"), ]
    if (dname == "OATS_MISS") d <- d[!(d$Block == "II" & d$Variety == "Victory" & d$Nitro == "0.6"), ]
    d <- droplevels(d)
    expect_equal(nrow(d), as.integer(S$est[S$table == "nobs"][1]))
    m <- lme4::lmer(Yield ~ Variety * Nitro + (1 | Block) + (1 | Block:Variety), data = d)
    dfm <- c(CONTAIN = "containment", BETWITHIN = "between-within", SATTERTHWAITE = "Satterthwaite",
             KR = "Kenward-Roger")[[ddfm]]
    vc <- as.data.frame(lme4::VarCorr(m))
    cp <- S[S$table == "covparms", ]
    rcp <- c(Intercept = vc$vcov[vc$grp == "Block"], Variety = vc$vcov[vc$grp == "Block:Variety"],
             Residual = vc$vcov[vc$grp == "Residual"])
    add("covparm", cp$effect, cp$est, rcp[cp$effect])
    terms <- c(Variety = "Variety", Nitro = "Nitro", "Variety*Nitro" = "Variety:Nitro")
    jt <- rbind(pm_joint(m, "Variety", df_method = dfm), pm_joint(m, "Nitro", df_method = dfm),
                pm_joint(m, "Variety:Nitro", df_method = dfm,
                         ctrmatrix = t(kronecker(stats::contr.helmert(3), stats::contr.helmert(4)))))
  } else {
    d <- sas_rep()
    m <- nlme::gls(Y ~ Group * Time, data = d, correlation = nlme::corSymm(form = ~ as.integer(Time) | Subj),
                   weights = nlme::varIdent(form = ~ 1 | Time))
    dfm <- c(CONTAIN = "auto", BETWITHIN = "between-within", SATTERTHWAITE = "Satterthwaite", KR = "Kenward-Roger")[[ddfm]]
    V <- nlme::getVarCov(m)
    cp <- S[S$table == "covparms", ]
    ij <- regmatches(cp$effect, regexpr("[0-9],[0-9]", cp$effect))
    add("covparm", cp$effect, cp$est,
        vapply(strsplit(ij, ","), function(k) V[as.integer(k[1]), as.integer(k[2])], numeric(1)))
    terms <- c(Group = "Group", Time = "Time", "Group*Time" = "Group:Time")
    jt <- rbind(pm_joint(m, "Group", df_method = dfm), pm_joint(m, "Time", df_method = dfm),
                pm_joint(m, "Group:Time", df_method = dfm,
                         ctrmatrix = t(kronecker(stats::contr.helmert(2), stats::contr.helmert(3)))))
  }
  for (eff in names(terms)) {
    mt <- terms[[eff]]
    two <- grepl(":", mt)
    pm <- predictmeans(m, mt, df_method = dfm, plot = FALSE, prtnum = FALSE)
    L <- S[S$table == "lsmeans" & S$effect == eff, ]
    key <- if (two) paste(L$level1, L$level2, sep = ":") else L$level1
    mm <- pm$means
    vars <- strsplit(mt, ":")[[1]]
    mkey <- do.call(paste, c(lapply(mm[vars], as.character), sep = ":"))
    i <- match(key, mkey)
    add("lsmean est", key, L$est, mm$estimate[i])
    add("lsmean se", key, L$se, mm$std.error[i])
    add("lsmean df", key, L$df, mm$df[i])
    D <- S[S$table == "diffs" & S$effect == eff, ]
    k <- pair_keys(D, two)
    pr <- pm$pairs
    j <- match(paste(k[[1]], k[[2]]), paste(pr$level1, pr$level2))
    sgn <- ifelse(is.na(j), -1, 1)
    j[is.na(j)] <- match(paste(k[[2]], k[[1]]), paste(pr$level1, pr$level2))[is.na(j)]
    dk <- paste(k[[1]], "-", k[[2]])
    add("diff est", dk, D$est, sgn * pr$estimate[j])
    add("diff se", dk, D$se, pr$std.error[j])
    add("diff df", dk, D$df, pr$df[j])
  }
  T3 <- S[S$table == "tests3", ]
  i <- match(T3$effect, names(terms))
  add("F", T3$effect, T3$F, jt$F.ratio[i])
  add("den df", T3$effect, T3$dendf, jt$df2[i])
  out <- do.call(rbind, res)
  out$run <- paste(dname, ddfm)
  out
}

test_that("all 16 PROC MIXED runs agree with SAS (containment, between-within, Satterthwaite, Kenward-Roger)", {
  skip_on_cran()
  ref <- sas_ref()
  all <- do.call(rbind, lapply(split(ref, as.integer(ref$run)), function(S) suppressMessages(suppressWarnings(sas_run(S)))))
  all$ok <- sas_agree(all$r, all$sas)
  # LS-mean SEs depend on where the REML iterations stop: 2e-5 between SAS and lme4 here
  se <- all$what == "lsmean se"
  all$ok[se] <- all$ok[se] | abs(all$r[se] - as.numeric(all$sas[se])) <= 1e-4 * as.numeric(all$sas[se])
  bad <- all[!all$ok, ]
  expect(nrow(bad) == 0, paste0(nrow(bad), " values differ from SAS:\n",
                                paste(utils::capture.output(print(utils::head(bad, 20), row.names = FALSE)),
                                      collapse = "\n")))
  expect_gt(nrow(all), 2500)
  # every kind of value was compared, for every run
  expect_setequal(unique(all$what), c("covparm", "lsmean est", "lsmean se", "lsmean df", "diff est", "diff se",
                                      "diff df", "F", "den df"))
  expect_length(unique(all$run), 16)
})

test_that("pbkrtest's Kenward-Roger (KR-expected) differs from SAS's when a whole plot is missing; Kenward-Roger does not", {
  skip_on_cran()
  ref <- sas_ref()
  d <- droplevels(sas_oats()[!(sas_oats()$Block == "I" & sas_oats()$Variety == "GoldenRain"), ])
  m <- lme4::lmer(Yield ~ Variety * Nitro + (1 | Block) + (1 | Block:Variety), data = d)
  D <- ref[ref$ddfm == "KR" & ref$data == "OATS_MISSWP" & ref$table == "diffs" & ref$effect == "Variety", ]
  s <- predictmeans(m, "Variety", df_method = "Kenward-Roger", plot = FALSE, prtnum = FALSE)$pairs$df
  k <- predictmeans(m, "Variety", df_method = "KR-expected", plot = FALSE, prtnum = FALSE)$pairs$df
  expect_equal(s, as.numeric(D$df), tolerance = 2e-3)                      # 8.35, 8.35, 8.1
  expect_true(all(k / as.numeric(D$df) > 1.08))                            # 9.25, 9.25, 9.01
})
