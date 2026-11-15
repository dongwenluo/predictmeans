## Comparison intervals: intervals around the means whose overlap answers the
## pairwise tests.
##
## Two means differ at the level of the test when |d| > crit x SED. Their
## intervals mean +/- c SE do not overlap when |d| > c (SE1 + SE2). The two
## agree for
##     c = crit x SED / (SE1 + SE2),   SED^2 = SE1^2 + SE2^2 - 2 cov(mean1, mean2).
## Independent means with equal SEs: c = crit/sqrt(2), the 83.4% interval for
## a 5% test (Goldstein and Healy 1995; Austin and Hux 2002; the 1.39 of the
## notched box plot, Velleman and Hoaglin 1981). Independent means with SE
## ratio K: c = crit sqrt(1 + K^2)/(1 + K), from 83.4% towards 95%
## (Maghsoodloo and Huang 2010). Correlated means (blocks, repeated measures,
## mixed models) have a smaller SED and so a smaller c (Afshartous and Preston
## 2010): the 83.4% rule is then too wide, and 95% confidence intervals far
## too wide, to be read as tests.
##
## One multiplier serves all the means: the average of c over the pairs
## compared (all pairs, or the pairs within a level of 'atvar'). When every
## pair has the same c (equal SEs and SEDs, as in balanced designs) overlap
## and test agree exactly; otherwise for most pairs, and 'disagree' counts
## those whose test says otherwise. 'crit' is the critical value of the
## test: t with the pair's df, the value that goes with 'adj' when it has one
## (tukey, bonferroni, sidak, scheffe, mvt), or the permutation quantile.
##
## pair_tab: pm_pair_table(); se, df: SEs and df of the means; level: the
## significance level. Returns NULL (no pair to compare) or list(mult,
## conf.level (the confidence level of mean +/- mult SE as an interval for
## one mean), level, adjust, adjusted, pairs, disagree, exact, range).
pm_comparison_intervals <- function(pair_tab, se, df, level, adj = "none") {
  if (is.null(pair_tab) || !nrow(pair_tab)) return(NULL)
  within <- if (any(!is.na(pair_tab$by))) !is.na(pair_tab$by) else rep(TRUE, nrow(pair_tab))
  adjusted <- !is.null(pair_tab$crit_adj)
  crit <- if (adjusted) pair_tab$crit_adj else pair_tab$crit
  sum_se <- se[pair_tab$i] + se[pair_tab$j]
  cij <- crit * pair_tab$ci_se / sum_se
  ok <- within & is.finite(cij) & !is.na(pair_tab$estimate) & sum_se > 0
  if (!any(ok)) return(NULL)
  mult <- mean(cij[ok])
  d <- abs(pair_tab$estimate[ok])
  by_test <- d > crit[ok] * pair_tab$ci_se[ok]
  by_overlap <- d > mult * sum_se[ok]
  dfm <- if (is.null(df) || all(is.na(df))) Inf else mean(df, na.rm = TRUE)
  rg <- range(cij[ok])
  list(mult = mult, conf.level = 1 - 2 * stats::pt(-mult, dfm), level = level, adjust = adj, adjusted = adjusted,
       pairs = sum(ok), disagree = sum(by_test != by_overlap), exact = diff(rg) <= 1e-6 * mult, range = rg)
}

## What to say when overlap and tests do not agree for every pair.
pm_comparison_caveat <- function(cmp) {
  if (isTRUE(cmp$exact)) return("")
  if (cmp$disagree > 0) {
    return(paste0(" (the SEDs differ: ", cmp$disagree, " of ", cmp$pairs, " pairs are judged differently by their own test)"))
  }
  " (the SEDs differ: approximate, but every pair here agrees with its test)"
}
