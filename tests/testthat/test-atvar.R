oats_fit <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  d
}

test_that("atvar is checked early: misspelt names, all factors, 'A:B' and NULL spellings", {
  d <- oats_fit()
  fm <- lm(yield ~ nitro * Variety, data = d)
  expect_error(predictmeans(fm, "nitro:Variety", atvar = "Varety"), "must name factors of modelterm")
  expect_error(predictmeans(fm, "nitro:Variety", atvar = c("nitro", "Variety")), "no pairs would be left")
  a <- predictmeans(fm, "nitro:Variety", atvar = "Variety")
  expect_identical(a$pairs, predictmeans(fm, "nitro:Variety", atvar = c("Variety", "NULL"))$pairs)
  expect_null(predictmeans(fm, "nitro:Variety", atvar = "")$settings$atvar)
  expect_null(pm_check_atvar("Variety", "nitro"))
  expect_identical(pm_check_atvar("B:A", c("A", "B", "C")), c("B", "A"))
})

test_that("with atvar the LSDs are those of the pairs within its levels, with each family's critical value", {
  d <- oats_fit()
  set.seed(5172)
  d <- d[-sample(nrow(d), 9), ]                    # unequal SEDs
  fm <- lm(yield ~ nitro * Variety, data = d)
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
  fl <- attr(pm$LSD, "Within levels of atvar")
  w <- pm$pairs[!is.na(pm$pairs$by), ]
  lsd <- (w$adj.conf.high - w$adj.conf.low) / 2
  expect_equal(unname(pm$LSD[c("Max.LSD", "Min.LSD", "Aveg.LSD")]), c(max(lsd), min(lsd), mean(lsd)),
               tolerance = 1e-4)
  # Tukey HSD of 4 means with the residual df
  expect_equal(fl$crit, rep(qtukey(0.95, 4, df.residual(fm)) / sqrt(2), 3), tolerance = 1e-8)
  expect_equal(fl$pairs, rep(6L, 3))
  expect_equal(fl$max.LSD[1], max(lsd[w$by == fl$by[1]]), tolerance = 1e-8)
  # pairs across levels are not in any LSD
  expect_true(all(is.na(pm$pairs$by[-seq_len(18)])))
})

test_that("means and pairs follow atvar's order; meandecr sorts within its levels", {
  d <- oats_fit()
  fm <- lm(yield ~ nitro * Variety, data = d)
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", meandecr = TRUE)
  expect_identical(as.character(pm$means$Variety), rep(levels(d$Variety), each = 4))
  expect_identical(unique(pm$pairs$by[!is.na(pm$pairs$by)]), levels(d$Variety))
  expect_equal(unname(coef(pm)), pm$means$estimate)
  tab <- pm$mean_table$Table
  expect_identical(names(tab)[1:2], c("Variety", "nitro"))
  for (v in levels(d$Variety)) expect_false(is.unsorted(rev(tab$Mean[tab$Variety == v])))
})

test_that("back-transformed columns stay with their means when the table is reordered", {
  d <- oats_fit()
  fm <- lm(log(yield) ~ nitro * Variety, data = d)
  pm <- suppressMessages(predictmeans(fm, "nitro:Variety", atvar = "Variety", meandecr = TRUE))
  tab <- pm$mean_table$Table
  expect_equal(tab$Bk_Mean, round(exp(tab$Mean), 4), tolerance = 1e-3)
})
