oats_joint <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  set.seed(7)
  d[-sample(nrow(d), 9), ]
}

test_that("pm_joint() slices: Kenward-Roger equals pbkrtest, Satterthwaite equals lmerTest", {
  skip_if_not_installed("pbkrtest")
  d <- oats_joint()
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = d)
  km <- Kmatrix(fm, "nitro:Variety")
  C <- contrast_family(km, "consec", by = "Variety")
  j <- pm_joint(fm, "nitro:Variety", by = "Variety", df_method = "KR-expected")
  js <- pm_joint(fm, "nitro:Variety", by = "Variety", df_method = "Satterthwaite")
  expect_equal(as.character(j$Variety), levels(d$Variety))
  for (g in seq_len(3)) {
    L <- (C %*% km$K)[grepl(paste0("^", j$Variety[g], ": "), rownames(C)), ]
    kr <- pbkrtest::KRmodcomp(fm, L)$test["Ftest", ]
    expect_equal(c(j$df1[g], j$df2[g], j$F.ratio[g], j$p.value[g]),
                 c(kr$ndf, kr$ddf, kr$stat, kr$p.value), tolerance = 1e-8)
    st <- lmerTest::contest(lmerTest::as_lmerModLmerTest(fm), L, joint = TRUE, ddf = "Satterthwaite")
    expect_equal(c(js$df2[g], js$F.ratio[g], js$p.value[g]), c(st$DenDF, st$`F value`, st$`Pr(>F)`), tolerance = 1e-6)
  }
  expect_identical(attr(j, "df_method"), "KR-expected")
})

test_that("pm_joint() for lm equals the F test of dropping the term; dependent rows are dropped", {
  d <- oats_joint()
  fl <- lm(yield ~ Block + nitro + Variety, data = d)
  j <- pm_joint(fl, "nitro")
  a <- anova(update(fl, . ~ . - nitro), fl)
  expect_equal(c(j$df1, j$df2, j$F.ratio, j$p.value), c(a$Df[2], a$Res.Df[2], a$F[2], a$`Pr(>F)`[2]), tolerance = 1e-10)
  cm <- contrast_family(Kmatrix(fl, "nitro"), "pairwise")       # 6 pairs, rank 3
  j2 <- pm_joint(fl, "nitro", ctrmatrix = cm)
  expect_equal(j2$df1, 3)
  expect_equal(j2$F.ratio, j$F.ratio, tolerance = 1e-10)
})

test_that("pm_joint() for lme: Kenward-Roger by default; containment is nlme's df of the term, as anova()", {
  d <- oats_joint()
  fm <- nlme::lme(yield ~ nitro*Variety, random = ~1|Block/Variety, data = d)
  j <- pm_joint(fm, "nitro:Variety", by = "nitro", df_method = "containment")
  expect_equal(unique(j$df2), unname(terms(fm$fixDF)["nitro:Variety"]))
  expect_equal(nrow(j), 4)
  expect_identical(attr(j, "df_method"), "containment (nlme)")
  expect_equal(pm_joint(fm, "Variety", df_method = "containment")$df2, unname(terms(fm$fixDF)["Variety"]))
  pm <- predictmeans(fm, "nitro:Variety", df_method = "containment", plot = FALSE, prtnum = FALSE)
  expect_equal(unique(pm$means$df), unique(j$df2))
  expect_identical(pm$settings$df_method, attr(j, "df_method"))
  expect_identical(attr(pm_joint(fm, "Variety"), "df_method"), "Kenward-Roger")
})
