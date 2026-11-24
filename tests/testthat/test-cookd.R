test_that("CookD(): lm and aov give the exact closed form; deletion refits agree with it", {
  O <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  fl <- lm(yield ~ Block + nitro*Variety, O)
  c1 <- CookD(fl, plot = FALSE)
  expect_equal(unname(c1), unname(cooks.distance(fl)), tolerance = 1e-10)
  expect_equal(CookD(aov(yield ~ Block + nitro*Variety, O), plot = FALSE), c1)
  # the definition: (b - b(-i))' V^-1 (b - b(-i)) / p, with V of the full fit
  b <- coef(fl); V <- vcov(fl)
  for (i in c(3, 40)) {
    bi <- coef(lm(yield ~ Block + nitro*Variety, O[-i, ]))
    expect_equal(unname(c1[i]), drop(t(b - bi) %*% solve(V) %*% (b - bi))/length(b), tolerance = 1e-8)
  }
  # deletion by group: deleting a whole block leaves its fixed effect without data: NA
  O$Block <- factor(O$Block, ordered = FALSE)
  fb <- lm(yield ~ Block + nitro*Variety, O)
  cg <- CookD(fb, group = "Block", plot = FALSE)
  expect_length(cg, 6)
  expect_true(is.na(cg[["II"]]))
  # by variety within a model without variety effects: the definition
  fv <- lm(yield ~ Block + nitro, O)
  cv <- CookD(fv, group = "Variety", plot = FALSE)
  b <- coef(fv); bi <- coef(lm(yield ~ Block + nitro, O[O$Variety != "Victory", ]))
  expect_equal(cv[["Victory"]], drop(t(b - bi) %*% solve(vcov(fv)) %*% (b - bi))/length(b), tolerance = 1e-8)
})

test_that("CookD(): glm by exact refits (cooks.distance() for a glm is a one-step approximation)", {
  g <- glm(breaks ~ wool*tension, poisson, warpbreaks)
  cd <- CookD(g, plot = FALSE)
  expect_length(cd, nrow(warpbreaks))
  b <- coef(g); V <- vcov(g)
  i <- which.max(cd)
  bi <- coef(glm(breaks ~ wool*tension, poisson, warpbreaks[-i, ]))
  expect_equal(unname(cd[i]), drop(t(b - bi) %*% solve(V) %*% (b - bi))/length(b), tolerance = 1e-8)
  # close to the one-step approximation, not equal
  expect_equal(unname(cd), unname(cooks.distance(g)), tolerance = 0.2)
})

test_that("CookD(): glmer and glmmTMB by groups, as the same model by the other package", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  d <- lme4::cbpp
  g1 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = d, family = binomial)
  g2 <- glmmTMB::glmmTMB(cbind(incidence, size - incidence) ~ period + (1 | herd), data = d, family = binomial)
  c1 <- CookD(g1, group = "herd", plot = FALSE)
  c2 <- CookD(g2, group = "herd", plot = FALSE)
  expect_length(c1, nlevels(d$herd))
  expect_equal(c1, c2, tolerance = 0.05)
  expect_equal(names(which.max(c1)), names(which.max(c2)))
})

test_that("CookD() refuses aov with Error() and other classes, with a clear message", {
  O <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  expect_error(CookD(aov(yield ~ nitro + Error(Block), O)), "Error\\(\\)")
  expect_error(CookD(list(a = 1)), "handles lm, aov, glm")
})

test_that("CookD() draws a ggplot with the largest distances labelled", {
  O <- as.data.frame(nlme::Oats); O$nitro <- factor(O$nitro)
  fl <- lm(yield ~ Block + nitro * Variety, O)
  cd <- CookD(fl, plot = FALSE)
  p <- pm_cookd_plot(cd, id = 2)
  expect_s3_class(p, "ggplot")
  lab <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)[[1]]$data
  expect_setequal(lab$label, names(sort(cd, decreasing = TRUE))[1:2])
  expect_no_error(ggplot2::ggplot_build(p))
  # group deletions: the groups on the x axis
  pg <- pm_cookd_plot(CookD(fl, group = "Block", plot = FALSE), group = "Block")
  expect_identical(pg$labels$x, "Level of Block")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_equal(CookD(fl), cd)                                   # drawing does not change the result
})
