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
## A mean over empty cells (R/empty-cells.R) has none when one of its pairs
## is not the difference of the shown means: compared on the cells both means
## have, or not estimable (no cell in common; confounded factors). No
## interval around such a mean can agree with its tests. Its pairs are left
## out, and the mean is named in 'left_out'.
##
## pair_tab: pm_pair_table(); est, se, df: the means, their SEs and df;
## level: the significance level; affected: TRUE for the means over empty
## cells. Returns NULL (no pair to compare) or list(mult, conf.level (the
## confidence level of mean +/- mult SE as an interval for one mean), level,
## adjust, adjusted, pairs, disagree (pairs whose overlap of the intervals
## around the shown means and test disagree), exact, range, left_out).
pm_comparison_intervals <- function(pair_tab, est, se, df, level, adj = "none", affected = NULL) {
  if (is.null(pair_tab) || !nrow(pair_tab)) return(NULL)
  if (is.null(affected)) affected <- rep(FALSE, length(se))
  within <- if (any(!is.na(pair_tab$by))) !is.na(pair_tab$by) else rep(TRUE, nrow(pair_tab))
  adjusted <- !is.null(pair_tab$crit_adj)
  crit <- if (adjusted) pair_tab$crit_adj else pair_tab$crit
  sum_se <- se[pair_tab$i] + se[pair_tab$j]
  cij <- crit * pair_tab$ci_se / sum_se
  off <- pm_misread_means(pair_tab, est, affected)
  ok <- within & is.finite(cij) & !is.na(pair_tab$estimate) & sum_se > 0 &
    !off[pair_tab$i] & !off[pair_tab$j] & is.finite(est[pair_tab$i]) & is.finite(est[pair_tab$j])
  if (!any(ok)) return(NULL)
  mult <- mean(cij[ok])
  by_test <- abs(pair_tab$estimate[ok]) > crit[ok] * pair_tab$ci_se[ok]
  # overlap is read from the means shown
  by_overlap <- abs(est[pair_tab$i[ok]] - est[pair_tab$j[ok]]) > mult * sum_se[ok]
  dfm <- if (is.null(df) || all(is.na(df))) Inf else mean(df, na.rm = TRUE)
  rg <- range(cij[ok])
  list(mult = mult, conf.level = 1 - 2 * stats::pt(-mult, dfm), level = level, adjust = adj, adjusted = adjusted,
       pairs = sum(ok), disagree = sum(by_test != by_overlap), exact = diff(rg) <= 1e-6 * mult, range = rg,
       left_out = names(affected)[off])
}

## The means over empty cells ('affected') whose pairs cannot be read off the
## shown means: TRUE for a mean with a pair (among those compared) that is not
## the difference of the shown means (compared on common cells) or is NA (no
## cell in common; confounded factors). The one rule for the comparison
## intervals, letterCI and the LSD note (pm_compare()).
pm_misread_means <- function(pair_tab, est, affected) {
  if (is.null(pair_tab) || !nrow(pair_tab) || !any(affected)) return(affected & FALSE)
  within <- if (any(!is.na(pair_tab$by))) !is.na(pair_tab$by) else rep(TRUE, nrow(pair_tab))
  shown <- est[pair_tab$j] - est[pair_tab$i]
  tol <- 1e-8 * max(1, abs(est), na.rm = TRUE)
  differs <- within & (affected[pair_tab$i] | affected[pair_tab$j]) &
    (is.na(pair_tab$estimate) | !is.finite(shown) |
       pmin(abs(pair_tab$estimate - shown), abs(pair_tab$estimate + shown)) > tol)
  affected & seq_along(affected) %in% c(pair_tab$i[differs], pair_tab$j[differs])
}

## What to say when overlap and tests do not agree for every pair.
pm_comparison_caveat <- function(cmp) {
  lo <- if (length(cmp$left_out)) paste0("; none for ", paste(cmp$left_out, collapse = ", "),
                                         " (means over empty cells: see the pairs)")
  if (isTRUE(cmp$exact)) return(if (length(lo)) paste0(" (", substring(lo, 3), ")") else "")
  if (cmp$disagree > 0) {
    return(paste0(" (the SEDs differ: ", cmp$disagree, " of ", cmp$pairs, " pairs are judged differently by their own test",
                  lo, ")"))
  }
  paste0(" (the SEDs differ: approximate, but every pair here agrees with its test", lo, ")")
}
