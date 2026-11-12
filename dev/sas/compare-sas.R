## Compare predictmeans with the SAS PROC MIXED reference output
## (dev/sas/sas-results.csv, parsed from the SAS Studio results of
## pm-reference.sas by parse-sas-results.py).
##
## Run: Rscript dev/sas/compare-sas.R [library] [csv]   (csv: default dev/sas/sas-results.csv)
## Prints, for each run and table, how many values agree with SAS to the
## precision SAS printed (half a unit in the last printed digit), and lists
## the values that do not.

lib <- commandArgs(TRUE)[1]
if (!is.na(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({ library(predictmeans); library(lme4); library(nlme) })
csv_file <- commandArgs(TRUE)[2]
if (is.na(csv_file)) csv_file <- "dev/sas/sas-results.csv"
sas <- read.csv(csv_file, colClasses = "character")

## data, as in pm-reference.sas (SAS reads Variety as 8 characters)
tok <- scan("dev/sas/oats.txt", what = "", quiet = TRUE)
oats <- data.frame(Block = tok[c(TRUE, FALSE, FALSE, FALSE)], Variety = substr(tok[c(FALSE, TRUE, FALSE, FALSE)], 1, 8),
                   Nitro = tok[c(FALSE, FALSE, TRUE, FALSE)], Yield = as.numeric(tok[c(FALSE, FALSE, FALSE, TRUE)]))
oats[1:3] <- lapply(oats[1:3], factor)
## The SAS run of 2026-11-12 deleted nothing for OATS_MISSWP: the script
## compared Variety (read as 8 characters, "GoldenRa") with 'GoldenRain'. So
## OATS_MISSWP was the full data (72 observations, as its "Number of
## Observations Used" shows) and OATS_MISS lacked only Block II x Victory x
## 0.6 (71). The data below are those SAS used; pm-reference.sas is fixed
## for later runs.
data_sets <- list(
  OATS = oats,
  OATS_MISSWP = oats,
  OATS_MISS = droplevels(subset(oats, !(Block == "II" & Variety == "Victory" & Nitro == "0.6"))))
rp <- read.table("dev/sas/rep.txt", col.names = c("Subj", "Group", "Time1", "Time2", "Time3"))
rep_long <- data.frame(Subj = factor(rep(rp$Subj, each = 3)), Group = factor(rep(rp$Group, each = 3)),
                       Time = factor(rep(1:3, nrow(rp))), Y = as.vector(t(as.matrix(rp[3:5]))))
rep_long <- rep_long[order(rep_long$Subj, rep_long$Time), ]

## agreement with the value SAS printed: within one unit of its last digit
## (half a unit is rounding; the rest allows for where the REML iterations
## stopped, a relative difference of about 1e-5)
decimals <- function(s) ifelse(grepl("\\.", s), nchar(sub("^[^.]*\\.", "", s)), 0)
agree <- function(r, s) {
  ok <- s != "" & !grepl("<", s)
  out <- rep(NA, length(s))
  out[ok] <- abs(r[ok] - as.numeric(s[ok])) <= 10^(-decimals(s[ok])) * (1 + 1e-6) + 1e-9
  out
}
report <- list()
check <- function(run, table, what, r, s, keys) {
  a <- agree(r, s)
  bad <- which(!is.na(a) & !a)
  report[[length(report) + 1]] <<- data.frame(run = run, table = table, what = what, n = sum(!is.na(a)),
                                              agree = sum(a, na.rm = TRUE))
  if (length(bad)) {
    cat(sprintf("  run %2s %-8s %-6s %-28s SAS %-10s R %s\n", run, table, what, keys[bad], s[bad],
                formatC(r[bad], digits = 6, format = "g")), sep = "")
  }
}
ddfm_r <- c(CONTAIN = "containment", BETWITHIN = "between-within", SATTERTHWAITE = "Satterthwaite",
            KR = "KR-SAS")       # "Kenward-Roger" (pbkrtest's expected information) differs by ~1% when unbalanced
lab <- function(x) substr(sub("^\\s+|\\s+$", "", x), 1, 8)   # Variety: first 8 characters, as in the first run

for (run in unique(sas$run)) {
  S <- sas[sas$run == run, ]
  dname <- S$data[1]
  ddfm <- S$ddfm[1]
  cat(sprintf("run %s: %s, DDFM=%s\n", run, dname, ddfm))
  if (dname != "REPEATEDLONG") {
    d <- data_sets[[dname]]
    m <- lmer(Yield ~ Variety * Nitro + (1 | Block) + (1 | Block:Variety), data = d, REML = TRUE)
    vc <- as.data.frame(VarCorr(m))
    cp <- S[S$table == "covparms", ]
    r_cp <- c(Intercept = vc$vcov[vc$grp == "Block"], Variety = vc$vcov[vc$grp == "Block:Variety"],
              Residual = vc$vcov[vc$grp == "Residual"])
    check(run, "covparm", "est", r_cp[cp$effect], cp$est, cp$effect)
    for (term in c("Variety", "Nitro", "Variety*Nitro")) {
      mt <- sub("\\*", ":", term)
      pm <- suppressMessages(suppressWarnings(
        predictmeans(m, mt, df_method = ddfm_r[[ddfm]], plot = FALSE, prtnum = FALSE, pairwise = TRUE)))
      L <- S[S$table == "lsmeans" & S$effect == term, ]
      key <- if (term == "Variety*Nitro") paste(lab(L$level1), lab(L$level2), sep = ":") else lab(L$level1)
      mm <- pm$means
      mkey <- if (term == "Variety*Nitro") paste(mm$Variety, mm$Nitro, sep = ":") else as.character(mm[[mt]])
      i <- match(key, mkey)
      check(run, "lsmeans", "est", mm$estimate[i], L$est, key)
      check(run, "lsmeans", "se", mm$std.error[i], L$se, key)
      check(run, "lsmeans", "df", mm$df[i], L$df, key)
      D <- S[S$table == "diffs" & S$effect == term, ]
      k1 <- if (term == "Variety*Nitro") paste(lab(D$level1), lab(D$level2), sep = ":") else lab(D$level1)
      k2 <- if (term == "Variety*Nitro") paste(lab(D$level3), lab(D$level4), sep = ":") else lab(D$level2)
      pr <- pm$pairs
      j <- match(paste(k1, k2), paste(pr$level1, pr$level2))
      sgn <- rep(1, length(j))
      if (anyNA(j)) {          # the pair in the other order
        j2 <- match(paste(k2, k1), paste(pr$level1, pr$level2))
        sgn[is.na(j)] <- -1
        j[is.na(j)] <- j2[is.na(j)]
      }
      dk <- paste(k1, "-", k2)
      check(run, "diffs", "est", sgn * pr$estimate[j], D$est, dk)
      check(run, "diffs", "se", pr$std.error[j], D$se, dk)
      check(run, "diffs", "df", pr$df[j], D$df, dk)
    }
    # Type 3 tests: equality of the LS-means (main effects), interaction contrasts
    T3 <- S[S$table == "tests3", ]
    cV <- contr.helmert(nlevels(d$Variety)); cN <- contr.helmert(nlevels(d$Nitro))
    jt <- rbind(
      suppressWarnings(pm_joint(m, "Variety", df_method = ddfm_r[[ddfm]])),
      suppressWarnings(pm_joint(m, "Nitro", df_method = ddfm_r[[ddfm]])),
      suppressWarnings(pm_joint(m, "Variety:Nitro", ctrmatrix = t(kronecker(cV, cN)), df_method = ddfm_r[[ddfm]])))
    rownames(jt) <- c("Variety", "Nitro", "Variety*Nitro")
    i <- match(T3$effect, rownames(jt))
    check(run, "tests3", "F", jt$F.ratio[i], T3$F, T3$effect)
    check(run, "tests3", "dendf", jt$df2[i], T3$dendf, T3$effect)
  } else {
    g <- gls(Y ~ Group * Time, data = rep_long, method = "REML",
             correlation = corSymm(form = ~ as.integer(Time) | Subj), weights = varIdent(form = ~ 1 | Time))
    V <- getVarCov(g)
    cp <- S[S$table == "covparms", ]
    ij <- regmatches(cp$effect, regexpr("[0-9],[0-9]", cp$effect))
    r_cp <- vapply(strsplit(ij, ","), function(k) V[as.integer(k[1]), as.integer(k[2])], numeric(1))
    check(run, "covparm", "est", r_cp, cp$est, cp$effect)
    for (term in c("Group", "Time", "Group*Time")) {
      mt <- sub("\\*", ":", term)
      dfm <- if (ddfm == "BETWITHIN") "between-within" else "auto"
      pm <- suppressMessages(suppressWarnings(predictmeans(g, mt, plot = FALSE, prtnum = FALSE, pairwise = TRUE,
                                                           df_method = dfm)))
      L <- S[S$table == "lsmeans" & S$effect == term, ]
      key <- if (term == "Group*Time") paste(lab(L$level1), lab(L$level2), sep = ":") else lab(L$level1)
      mm <- pm$means
      mkey <- if (term == "Group*Time") paste(mm$Group, mm$Time, sep = ":") else as.character(mm[[mt]])
      i <- match(key, mkey)
      check(run, "lsmeans", "est", mm$estimate[i], L$est, key)
      check(run, "lsmeans", "se", mm$std.error[i], L$se, key)
      if (ddfm %in% c("CONTAIN", "BETWITHIN")) check(run, "lsmeans", "df", mm$df[i], L$df, key)
      D <- S[S$table == "diffs" & S$effect == term, ]
      k1 <- if (term == "Group*Time") paste(lab(D$level1), lab(D$level2), sep = ":") else lab(D$level1)
      k2 <- if (term == "Group*Time") paste(lab(D$level3), lab(D$level4), sep = ":") else lab(D$level2)
      pr <- pm$pairs
      j <- match(paste(k1, k2), paste(pr$level1, pr$level2))
      sgn <- ifelse(is.na(j), -1, 1)
      j[is.na(j)] <- match(paste(k2, k1), paste(pr$level1, pr$level2))[is.na(j)]
      dk <- paste(k1, "-", k2)
      check(run, "diffs", "est", sgn * pr$estimate[j], D$est, dk)
      check(run, "diffs", "se", pr$std.error[j], D$se, dk)
      if (ddfm %in% c("CONTAIN", "BETWITHIN")) check(run, "diffs", "df", pr$df[j], D$df, dk)
    }
  }
}
rep_df <- do.call(rbind, report)
agg <- aggregate(cbind(n, agree) ~ table + what, data = rep_df, FUN = sum)
cat("\nAgreement with SAS at its printed precision (all runs):\n")
print(agg[order(agg$table, agg$what), ], row.names = FALSE)
