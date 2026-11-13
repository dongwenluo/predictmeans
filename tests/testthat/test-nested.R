## Nested treatment structures (as in "Experimental Design and Analysis for
## Tree Improvement", Exams 4.3 and 8.1): combinations that do not exist by
## design are not empty cells. predictmeans averages over the combinations
## that exist, with equal weights, so in balanced data a mean equals the raw
## mean. emmeans weights by the nesting hierarchy instead (for Exam 4.3 it
## gives the single control treatment half the weight of a seedlot mean), so
## it is not the reference here.

test_that("treatments nested in a control/treated factor: means over the treatments that exist", {
  set.seed(43)
  d <- expand.grid(repl = factor(1:3), seedlot = factor(1:6), treat = factor(c("control", "nick", "bws", "bw1min"),
                                                                           levels = c("control", "nick", "bws", "bw1min")))
  d$contcomp <- factor(ifelse(d$treat == "control", "control", "Treated"))
  d$y <- 50 + 20 * (d$treat != "control") + as.integer(d$seedlot) + rnorm(nrow(d), 0, 3)
  fm <- aov(y ~ repl + contcomp + seedlot + treat/contcomp + contcomp/seedlot + treat/contcomp/seedlot, data = d)
  pm <- function(term) suppressWarnings(predictmeans(fm, term, plot = FALSE, prtnum = FALSE))$means
  s <- pm("seedlot")
  expect_true(all(s$estimable))
  expect_equal(s$estimate, as.vector(tapply(d$y, d$seedlot, mean)), tolerance = 1e-8)   # each treatment weight 1/4
  expect_equal(pm("treat")$estimate, as.vector(tapply(d$y, d$treat, mean)), tolerance = 1e-8)
  expect_equal(pm("contcomp")$estimate, as.vector(tapply(d$y, d$contcomp, mean)), tolerance = 1e-8)
  cs <- pm("contcomp:seedlot")
  expect_equal(nrow(cs), 12)
  expect_true(all(cs$estimable))
})

test_that("provenances nested in countries, crossed with a treatment: equal weights over provenances", {
  set.seed(81)
  pv <- data.frame(country = rep(c("A", "B", "C"), c(2, 3, 5)), prov = paste0("p", 1:10))
  d <- merge(expand.grid(prov = pv$prov, inoc = c("no", "yes"), repl = factor(1:4)), pv)
  d[] <- lapply(d, function(v) if (is.character(v)) factor(v) else v)
  d$dbh <- 10 + as.integer(d$country) + rnorm(10)[as.integer(d$prov)] + 0.5 * (d$inoc == "yes") + rnorm(nrow(d), 0, 0.5)
  fm <- lm(dbh ~ repl + inoc * (country/prov), data = d)
  pm <- function(term) suppressWarnings(predictmeans(fm, term, plot = FALSE, prtnum = FALSE))$means
  expect_equal(pm("inoc")$estimate, as.vector(tapply(d$dbh, d$inoc, mean)), tolerance = 1e-8)
  expect_equal(pm("country")$estimate, as.vector(tapply(d$dbh, d$country, mean)), tolerance = 1e-8)
  ic <- pm("inoc:country")
  expect_equal(ic$estimate, as.vector(t(tapply(d$dbh, list(d$inoc, d$country), mean))), tolerance = 1e-8)
})

test_that("an empty cell (a combination missing from the data): the 1.x average of the observed cells, with a warning", {
  set.seed(81)
  pv <- data.frame(country = rep(c("A", "B", "C"), c(2, 3, 5)), prov = paste0("p", 1:10))
  d <- merge(expand.grid(prov = pv$prov, inoc = c("no", "yes"), repl = factor(1:4)), pv)
  d[] <- lapply(d, function(v) if (is.character(v)) factor(v) else v)
  d$dbh <- 10 + as.integer(d$country) + rnorm(10)[as.integer(d$prov)] + 0.5 * (d$inoc == "yes") + rnorm(nrow(d), 0, 0.5)
  d <- d[!(d$prov == "p3" & d$inoc == "yes"), ]                      # p3 never received inoculation
  fm <- lm(dbh ~ repl + inoc * (country/prov), data = d)
  expect_warning(p <- predictmeans(fm, "inoc", plot = FALSE, prtnum = FALSE), "Missing treatments' combination")
  cells <- tapply(d$dbh, list(d$inoc, d$prov), mean)                  # repl is balanced within each cell
  expect_equal(p$means$estimate, unname(rowMeans(cells, na.rm = TRUE)), tolerance = 1e-8)
  expect_true(all(p$means$estimable))
})
