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
