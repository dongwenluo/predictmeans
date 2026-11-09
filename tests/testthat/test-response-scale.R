resp_data <- function() {
  d <- data.frame(A = gl(4, 10, labels = c("a1", "a2", "a3", "a4")), B = gl(2, 5, 40))
  set.seed(3)
  d$y <- rpois(40, c(3, 5, 8, 4)[d$A])
  d$s <- rbinom(40, 10, c(.2, .4, .5, .3)[d$A])
  d$z <- exp(rnorm(40, 1 + 0.2*as.numeric(d$A)))
  d
}

test_that("response-scale SEs (delta method) and ratios match emmeans", {
  skip_if_not_installed("emmeans")
  d <- resp_data()
  cases <- list(
    list(fit = glm(y ~ A + B, family = poisson, data = d), trans = NULL, col = "ratio"),
    list(fit = glm(cbind(s, 10 - s) ~ A + B, family = binomial, data = d), trans = NULL, col = "odds.ratio"),
    list(fit = lm(log(z) ~ A + B, data = d), trans = exp, col = "ratio"))
  for (cs in cases) {
    pm <- if (is.null(cs$trans)) predictmeans(cs$fit, "A", pairwise = TRUE, plot = FALSE) else
      predictmeans(cs$fit, "A", pairwise = TRUE, trans = cs$trans, plot = FALSE)
    em <- emmeans::emmeans(cs$fit, ~A)
    e <- summary(em, type = "response")
    i <- match(as.character(pm$means$A), as.character(e$A))
    expect_equal(pm$means$response, e[[2]][i], tolerance = 1e-10)
    expect_equal(pm$means$response.std.error, e$SE[i], tolerance = 1e-10)
    ep <- summary(pairs(em), type = "response", infer = TRUE, adjust = "none")
    j <- match(paste(pm$pairs$level1, "/", pm$pairs$level2), as.character(ep$contrast))
    expect_equal(pm$pairs[[cs$col]], ep[[2]][j], tolerance = 1e-10)
    # limits from the difference's CI, with the pairs' df
    crit <- qt(0.975, pm$pairs$df)
    expect_equal(pm$pairs[[paste0(cs$col, ".conf.low")]], exp(pm$pairs$estimate - crit*pm$pairs$std.error))
  }
  # lm on the log scale: same df as emmeans, so the limits agree too
  pm <- predictmeans(cases[[3]]$fit, "A", pairwise = TRUE, trans = exp, plot = FALSE)
  ep <- summary(pairs(emmeans::emmeans(cases[[3]]$fit, ~A)), type = "response", infer = TRUE, adjust = "none")
  j <- match(paste(pm$pairs$level1, "/", pm$pairs$level2), as.character(ep$contrast))
  expect_equal(pm$pairs$ratio.conf.low, ep$lower.CL[j], tolerance = 1e-10)
  expect_equal(pm$pairs$ratio.conf.high, ep$upper.CL[j], tolerance = 1e-10)
})

test_that("no response columns or ratios on the identity scale", {
  d <- resp_data()
  pm <- predictmeans(lm(z ~ A + B, data = d), "A", pairwise = TRUE, plot = FALSE)
  expect_null(pm$means$response.std.error)
  expect_false(any(grepl("ratio", names(pm$pairs))))
})
