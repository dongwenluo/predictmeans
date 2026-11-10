test_that("factors not in modelterm (incl. 'by' factors) are averaged in both parts", {
  skip_on_cran()
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  W <- WarsawApts
  W$district <- factor(W$district)
  sm <- suppressWarnings(semireg(areaPerMzloty ~ construction.date*district,
                                 smoothZ = list(grp = smZ(construction.date, k = 15, by = district, group = TRUE)),
                                 data = W, prt = FALSE))
  overall <- semipred(sm, covariate = "construction.date", prt = FALSE)$pred_df
  by_district <- semipred(sm, "district", "construction.date", prt = FALSE)$pred_df
  avg <- aggregate(by_district$areaPerMzloty, list(cd = by_district$construction.date), mean)
  expect_equal(overall$construction.date, avg$cd)
  expect_equal(overall$areaPerMzloty, avg$x, tolerance = 1e-10)
})

test_that("a matrix covariateV gives one prediction per point for a 2-D smooth", {
  set.seed(3)
  s2 <- data.frame(u = runif(200), v = runif(200))
  s2$w <- sin(3*s2$u) + cos(3*s2$v) + rnorm(200, sd = 0.2)
  sm <- semireg(w ~ u + v, smoothZ = list(g = smZ(cbind(u, v), k = 6, type = "Ztps")),
                data = s2, prt = FALSE)
  pts <- cbind(u = c(0.2, 0.5, 0.8), v = c(0.3, 0.5, 0.7))
  p <- semipred(sm, covariate = c("u", "v"), covariateV = pts, prt = FALSE)$pred_df
  expect_equal(nrow(p), 3)
  expect_equal(as.matrix(p[c("u", "v")]), pts, ignore_attr = TRUE, tolerance = 1e-12)
  expect_true(all(p$LL < p$w & p$w < p$UL))
})

test_that("a numeric variable that is not in the prediction (here the grouping id) doesn't stop semipred()", {
  set.seed(7)
  d <- data.frame(u = runif(150, 2, 8), id = rep(1:15, each = 10))
  d$y <- sin(d$u) + rnorm(15, sd = 0.5)[d$id] + rnorm(150, sd = 0.3)
  centred <- suppressMessages(semireg(y ~ u + (1 | id), smoothZ = list(g = smZ(u)), data = d, prt = FALSE))
  plain <- suppressMessages(semireg(y ~ u + (1 | id), smoothZ = list(g = smZ(u)), data = d,
                                    ncenter = FALSE, prt = FALSE))
  expect_equal(semipred(centred, covariate = "u", prt = FALSE)$pred_df,
               semipred(plain, covariate = "u", prt = FALSE)$pred_df, tolerance = 1e-6)
})

test_that("ranova() on a semireg fit refits with the smooth bases (copper)", {
  skip_on_cran()
  data(copper)
  full <- suppressWarnings(semireg(log(grade) ~ xcoord + ycoord + zcoord,
    smoothZ = list(grp_z = smZ(zcoord), grp_geo = smZ(cbind(xcoord, ycoord), type = "Ztps")),
    data = copper, prt = FALSE, predict_info = FALSE))
  expect_s4_class(full$semer, "semerLmerTest")
  rt <- suppressWarnings(lmerTest::ranova(full$semer))
  expect_equal(nrow(rt), 3)
  # dropping grp_geo must equal fitting with the grp_z smooth only
  only_z <- suppressWarnings(semireg(log(grade) ~ xcoord + ycoord + zcoord,
    smoothZ = list(grp_z = smZ(zcoord)), data = copper, prt = FALSE, predict_info = FALSE))
  expect_equal(rt$logLik[rownames(rt) == "(1 | grp_geo)"], as.numeric(logLik(only_z$semer)),
               tolerance = 1e-6)
})

test_that("predictmeans() on a semireg fit: means with the smooth, as semipred(); pairs as semipred(contr =)", {
  skip_on_cran()
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  W <- WarsawApts
  W$district <- factor(W$district)
  sm <- suppressWarnings(semireg(areaPerMzloty ~ construction.date*district,
                                 smoothZ = list(grp = smZ(construction.date, k = 15, by = district, group = TRUE)),
                                 data = W, prt = FALSE))
  r <- predictmeans(sm, "district", covariate = c(construction.date = 1980), adj = "tukey")
  expect_s3_class(r, "pm_semimeans")
  sp <- semipred(sm, "district", "construction.date", covariateV = 1980, prt = FALSE)$pred_df
  expect_equal(r$means$estimate, sp$areaPerMzloty, tolerance = 1e-10)
  expect_equal(r$means$std.error, sp$SE, tolerance = 1e-10)
  sc <- semipred(sm, "district", "construction.date", contr = c(2, 1), covariateV = 1980, prt = FALSE)$pred_df
  expect_equal(r$pairs$estimate[1], -sc$areaPerMzloty, tolerance = 1e-10)
  expect_equal(r$pairs$std.error[1], sc$SE_contr, tolerance = 1e-10)
  expect_equal(unname(sqrt(diag(vcov(r)))), r$means$std.error)
  # default: covariates at their means
  r0 <- predictmeans(sm, "district")
  sp0 <- semipred(sm, "district", "construction.date", covariateV = mean(W$construction.date), prt = FALSE)$pred_df
  expect_equal(r0$means$estimate, sp0$areaPerMzloty, tolerance = 1e-8)
})

test_that("semipred(): two-sided limits; simultaneous bands are wider, repeat, and keep the user's RNG", {
  skip_on_cran()
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  sm <- suppressWarnings(semireg(areaPerMzloty ~ construction.date,
                                 smoothZ = list(grp = smZ(construction.date, k = 25)), data = WarsawApts, prt = FALSE))
  p <- semipred(sm, covariate = "construction.date", prt = FALSE)$pred_df
  expect_equal((p$UL - p$areaPerMzloty)/p$SE, rep(qnorm(0.975), nrow(p)), tolerance = 1e-10)
  set.seed(9); u1 <- runif(1); set.seed(9)
  s <- semipred(sm, covariate = "construction.date", prt = FALSE, simultaneous = TRUE)$pred_df
  expect_identical(runif(1), u1)
  crit <- unique(round((s$UL - s$areaPerMzloty)/s$SE, 8))
  expect_length(crit, 1)
  expect_gt(crit, qnorm(0.975))
  expect_lt(crit, qnorm(1 - 0.05/(2*nrow(s))))      # below the Bonferroni bound
  s2 <- semipred(sm, covariate = "construction.date", prt = FALSE, simultaneous = TRUE)$pred_df
  expect_identical(s$UL, s2$UL)
})
