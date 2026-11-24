# A covariate that is constant within the levels of a factor (site mean in a
# Finlay-Wilkinson regression) is aliased with that factor; the means must use
# each level's own value, not the overall mean (which made every mean NA).

exam52 <- function() {
  sites <- c("A", "B", "C", "D")
  lots <- c("l1", "l2", "l3", "l4", "l5")
  set.seed(4817)
  d <- expand.grid(seedlot = lots, site = sites)
  eff <- c(A = 100, B = 300, C = 500, D = 200)
  d$ht <- eff[d$site] * c(0.8, 1, 1.1, 1.2, 0.9)[as.integer(d$seedlot)] + rnorm(nrow(d), 0, 10)
  d$ht[c(3, 9, 17)] <- NA
  d$sitemean <- ave(d$ht, d$site, FUN = function(z) mean(z, na.rm = TRUE))
  d$site <- factor(d$site)
  d$seedlot <- factor(d$seedlot)
  d
}

test_that("site-level covariate aliased with site gives estimable seedlot means", {
  d <- exam52()
  m <- lm(ht ~ site + sitemean * seedlot, d)
  expect_true(anyNA(coef(m)))
  expect_message(res <- predictmeans(m, "seedlot", plot = FALSE, weights = "cells"), "constant within the levels of site")
  mt <- res$mean_table
  expect_false(anyNA(mt$Mean))
  expect_false(anyNA(mt$SE))
  dd <- d[!is.na(d$ht), ]
  expect_equal(mt$Mean, unname(c(tapply(fitted(m), dd$seedlot, mean))), tolerance = 1e-6)

  eq <- suppressMessages(predictmeans(m, "seedlot", plot = FALSE))$mean_table
  expect_false(anyNA(eq$Mean))
  # a seedlot observed at every site: equal and cell weights agree
  full <- names(which(table(dd$seedlot) == 4))
  expect_equal(eq$Mean[eq$seedlot %in% full], mt$Mean[mt$seedlot %in% full], tolerance = 1e-6)
})

test_that("seedlot means do not depend on the reference level", {
  d <- exam52()
  m1 <- lm(ht ~ site + sitemean * seedlot, d)
  d2 <- d
  d2$seedlot <- relevel(d2$seedlot, "l3")
  m2 <- lm(ht ~ site + sitemean * seedlot, d2)
  for (w in c("equal", "cells")) {
    a <- suppressMessages(predictmeans(m1, "seedlot", plot = FALSE, weights = w))$mean_table
    b <- suppressMessages(predictmeans(m2, "seedlot", plot = FALSE, weights = w))$mean_table
    b <- b[match(a$seedlot, b$seedlot), ]
    expect_equal(a$Mean, b$Mean, tolerance = 1e-8)
    expect_equal(a$SE, b$SE, tolerance = 1e-8)
  }
})

test_that("site means use each site's own covariate value", {
  d <- exam52()
  m <- lm(ht ~ site + sitemean * seedlot, d)
  mt <- suppressMessages(predictmeans(m, "site", plot = FALSE))$mean_table
  expect_false(anyNA(mt$Mean))
  nd <- expand.grid(seedlot = levels(d$seedlot), site = levels(d$site))
  nd$sitemean <- tapply(d$sitemean, d$site, `[`, 1)[as.character(nd$site)]
  ref <- tapply(suppressWarnings(predict(m, nd)), nd$site, mean)
  expect_equal(mt$Mean, unname(c(ref)), tolerance = 1e-6)
})

# Williams, Matheson & Harwood, Experimental Design and Analysis for Tree
# Improvement, Table 5.7 (4 sites, 27 complete seedlots of Table 5.2).
test_that("tree data (DataExam5.2) agrees with the book's Table 5.7", {
  data("DataExam5.2", package = "predictmeans", envir = environment())
  d <- DataExam5.2

  # full, unbalanced data: no NA, cell means = averages of fitted values
  m <- lm(ht ~ site + sitemean * seedlot, d)
  mt <- suppressMessages(predictmeans(m, "seedlot", plot = FALSE, weights = "cells"))$mean_table
  expect_false(anyNA(mt$Mean))
  dd <- d[!is.na(d$ht), ]
  expect_equal(mt$Mean, unname(c(tapply(fitted(m), dd$seedlot, mean))), tolerance = 1e-6)

  # the book's balanced subset
  s <- d[d$site %in% c("Ratchaburi", "Sai Thong", "Si Sa Ket", "Sakaerat"), ]
  full <- names(which(tapply(!is.na(s$ht), s$seedlot, all)))
  s <- droplevels(s[s$seedlot %in% full, ])
  s$sitemean <- ave(s$ht, s$site)
  expect_equal(nlevels(s$seedlot), 27)
  m4 <- lm(ht ~ site + sitemean * seedlot, s)
  mt4 <- suppressMessages(predictmeans(m4, "seedlot", plot = FALSE))$mean_table
  book <- c(`13877` = 365, `13866` = 353, `13861` = 627, `13854` = 628, `13684` = 660, `13683` = 770,
            `12013` = 616, `14485` = 123, `11935` = 273, `14537` = 780)
  expect_equal(round(mt4$Mean[match(names(book), mt4$seedlot)]), unname(book))
  expect_equal(mt4$SE[1], 43.82294, tolerance = 1e-6)
  expect_equal(mt4$Df[1], 52)
  # the book's regression coefficients: each seedlot's slope on the site mean
  slope <- vapply(c("11935", "12013", "13877"),
                  function(l) coef(lm(ht ~ sitemean, s[s$seedlot == l, ]))[[2]], numeric(1))
  expect_equal(round(slope, 2), c(0.53, 2.07, 0.79), ignore_attr = TRUE)
})
