## Equivalence tests (TOST): predictmeans(equivalence =), contrastmeans(equivalence =)

oats_lm <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  lm(yield ~ Block + nitro + Variety, data = d)
}

test_that("TOST p-values of pairs: by hand and as emmeans", {
  fm <- oats_lm()
  pm <- predictmeans(fm, "nitro", equivalence = 20, plot = FALSE, prtnum = FALSE)
  pr <- pm$pairs
  byhand <- pmax(pt((pr$estimate + 20) / pr$std.error, pr$df, lower.tail = FALSE),
                 pt((pr$estimate - 20) / pr$std.error, pr$df))
  expect_equal(pr$p.equiv, byhand)
  expect_equal(pm$settings$equivalence, c(-20, 20))
  # asymmetric margins
  pa <- predictmeans(fm, "nitro", equivalence = c(-40, 5), plot = FALSE, prtnum = FALSE)$pairs
  expect_equal(pa$p.equiv, pmax(pt((pa$estimate + 40) / pa$std.error, pa$df, lower.tail = FALSE),
                                pt((pa$estimate - 5) / pa$std.error, pa$df)))
  expect_null(predictmeans(fm, "nitro", plot = FALSE, prtnum = FALSE)$pairs$p.equiv)
  skip_if_not_installed("emmeans")
  e <- summary(emmeans::test(pairs(emmeans::emmeans(fm, ~ nitro), adjust = "none"), delta = 20,
                             side = "equivalence"))
  key <- paste0("nitro", pr$level1, " - nitro", pr$level2)
  expect_equal(pr$p.equiv, e$p.value[match(key, as.character(e$contrast))], tolerance = 1e-8)
})

test_that("contrastmeans(equivalence =) and wrong margins", {
  fm <- oats_lm()
  cn <- rbind(c(-1, 1, 0, 0), c(0, 0, -1, 1))
  ct <- contrastmeans(fm, "nitro", ctrmatrix = cn, equivalence = 15)[[1]]
  est <- ct[, "Estimate"]; se <- ct[, "Std. Error"]; df <- ct[, "df"]
  expect_equal(unname(ct[, "p.equiv"]),
               unname(round(pmax(pt((est + 15) / se, df, lower.tail = FALSE), pt((est - 15) / se, df)), 4)),
               tolerance = 1e-3)
  expect_error(predictmeans(fm, "nitro", equivalence = -1, plot = FALSE, prtnum = FALSE), "equivalence")
  expect_error(contrastmeans(fm, "nitro", ctrmatrix = cn, equivalence = c(5, -5)), "equivalence")
})
