## predictmeans against SAS PROC MIXED: Tukey-Kramer adjusted differences of
## LS-means with Kenward-Roger df (dev/sas/pm-reference-adjust.sas, run in SAS
## OnDemand on 2026-11-13; parsed results in _data/sas/proc-mixed-adjust.csv.gz).
## Oats without the Block I x Golden Rain whole plot (68 observations): the
## pairs have different df (Variety 8.35 and 8.1; Variety x Nitro 42 and 22.2).
## predictmeans adjusts each pair with its own df, as SAS ADJDFE=ROW; SAS's
## default ADJDFE=SOURCE uses the denominator df of the term's F test for all
## its pairs (42 for Variety x Nitro, also for the pairs across whole plots
## whose own df are 22.2).

adj_fit <- function() {
  o <- as.data.frame(nlme::Oats)
  o$Variety <- factor(gsub(" ", "", o$Variety))
  o$Nitro <- factor(sprintf("%.1f", o$nitro))
  d <- droplevels(o[!(o$Block == "I" & o$Variety == "GoldenRain"), ])
  lme4::lmer(yield ~ Variety * Nitro + (1 | Block/Variety), data = d)
}
## within one unit of SAS's last printed digit ("<.0001": below)
adj_agree <- function(r, s) {
  dec <- ifelse(grepl("\\.", s), nchar(sub("^[^.]*\\.", "", s)), 0)
  v <- as.numeric(sub("<", "", s))
  ifelse(grepl("^<", s), r < v, abs(r - v) <= 10^(-dec) * (1 + 1e-6) + 1e-9)
}

test_that("adj = \"tukey\" equals SAS Tukey-Kramer with ADJDFE=ROW: estimates, SEs, df, p-values, limits", {
  skip_on_cran()
  S <- read.csv(test_path("_data", "sas", "proc-mixed-adjust.csv.gz"), colClasses = "character")
  fm <- adj_fit()
  pp <- function(term) suppressWarnings(predictmeans(fm, term, pairwise = TRUE, adj = "tukey", plot = FALSE,
                                                     prtnum = FALSE))$pairs
  P <- rbind(cbind(effect = "Variety", pp("Variety")[, c("level1", "level2", "estimate", "std.error", "df", "p.value",
                                                        "p.adj", "conf.low", "conf.high")]),
             cbind(effect = "Variety*Nitro", pp("Variety:Nitro")[, c("level1", "level2", "estimate", "std.error", "df",
                                                                     "p.value", "p.adj", "conf.low", "conf.high")]))
  for (run in c("row", "source")) {
    R <- S[S$adjdfe == run, ]
    l1 <- ifelse(R$effect == "Variety", R$variety, paste(R$variety, R$nitro, sep = ":"))
    l2 <- ifelse(R$effect == "Variety", R$variety2, paste(R$variety2, R$nitro2, sep = ":"))
    Q <- P[match(paste(R$effect, l1, l2), paste(P$effect, P$level1, P$level2)), ]
    expect_false(anyNA(Q$estimate))
    expect_true(all(adj_agree(Q$estimate, R$estimate)), label = paste(run, "estimates"))
    expect_true(all(adj_agree(Q$std.error, R$se)), label = paste(run, "SEs"))
    expect_true(all(adj_agree(Q$df, R$df)), label = paste(run, "df"))
    expect_true(all(adj_agree(Q$p.value, R$p)), label = paste(run, "p-values"))
    # confidence limits: within 1e-4 relative (the REML estimates of SAS and lme4
    # differ by about 1e-5 where the iterations stop)
    expect_equal(Q$conf.low, as.numeric(R$lower), tolerance = 1e-4)
    expect_equal(Q$conf.high, as.numeric(R$upper), tolerance = 1e-4)
    k <- ifelse(R$effect == "Variety", 3, 12)
    if (run == "row") {
      expect_true(all(adj_agree(Q$p.adj, R$adj_p)), label = "Tukey-Kramer p, ADJDFE=ROW")
      dfa <- Q$df
    } else {
      # ADJDFE=SOURCE: the same formula with the F test's denominator df
      dfa <- ifelse(R$effect == "Variety", pm_joint(fm, "Variety")$df2, 42)
      expect_true(all(adj_agree(ptukey(sqrt(2) * abs(Q$estimate / Q$std.error), k, dfa, lower.tail = FALSE), R$adj_p)),
                  label = "Tukey-Kramer p, ADJDFE=SOURCE")
    }
    crit <- qtukey(0.95, k, dfa) / sqrt(2)
    expect_equal(Q$estimate - crit * Q$std.error, as.numeric(R$adj_lower), tolerance = 1e-4)
    expect_equal(Q$estimate + crit * Q$std.error, as.numeric(R$adj_upper), tolerance = 1e-4)
  }
})
