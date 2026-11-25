## autoplot() of covariatemeans() and semipred() results (R/autoplot-curves.R)

q <- function(e) suppressWarnings(suppressMessages(e))
built <- function(p) { ggplot2::ggplot_build(p); TRUE }

test_that("covariatemeans(): one to four factors, default and chosen layouts", {
  skip_on_cran()
  O <- as.data.frame(nlme::Oats)
  O$B <- factor(O$Block %in% c("I", "II", "III"), labels = c("late", "early"))
  # C and D independent of nitro (by row position they would be functions of
  # nitro, and the curves not estimable)
  set.seed(1)
  O$C <- factor(sample(1:2, nrow(O), TRUE))
  O$D <- factor(sample(1:2, nrow(O), TRUE))
  fm <- lm(yield ~ nitro * Variety + B + C + D, O)
  terms <- c("Variety", "Variety:B", "Variety:B:C", "Variety:B:C:D")
  for (t in terms) {
    r <- q(covariatemeans(fm, t, covariate = "nitro", at = c(0, 0.3, 0.6)))
    p <- ggplot2::autoplot(r)
    expect_s3_class(p, "ggplot")
    expect_true(built(p))
    expect_s3_class(r$plt, "ggplot")                       # $plt: the same plot, for 1.x scripts
  }
  r2 <- q(covariatemeans(fm, "Variety:B", covariate = "nitro", at = c(0, 0.3, 0.6)))
  expect_s3_class(ggplot2::autoplot(r2)$facet, "FacetWrap")
  expect_s3_class(ggplot2::autoplot(r2, panel = c("Variety", "B"))$facet, "FacetGrid")
  p <- ggplot2::autoplot(r2, group = "B", panel = "Variety", data = FALSE, interval = "none")
  expect_equal(nlevels(ggplot2::ggplot_build(p)$layout$layout$PANEL), 3)
  expect_equal(p$labels$colour, "B")
  p1 <- ggplot2::autoplot(r2, panel = character(0))
  expect_s3_class(p1$facet, "FacetNull")
  expect_error(ggplot2::autoplot(r2, group = "Z"), "among the factors")
  # plot_args (superseded) map to the same arguments
  r3 <- q(covariatemeans(fm, "Variety:B", covariate = "nitro", at = c(0, 0.3, 0.6),
                         plot_args = list(plotord = c("B", "Variety"), ci = FALSE)))
  expect_equal(r3$plt$labels$colour, "B")
  expect_false("GeomRibbon" %in% vapply(r3$plt$layers, function(l) class(l$geom)[1], ""))
  expect_error(q(covariatemeans(fm, "Variety:B", covariate = "nitro", plot_args = list(plotord = 1))), "length 2")
})

test_that("covariatemeans(): pltdf is pred_df (message once), as_is keeps the right rows, response scale", {
  O <- as.data.frame(nlme::Oats)
  fm <- lm(log(yield) ~ nitro * Variety, O)
  r <- q(covariatemeans(fm, "Variety", covariate = "nitro"))
  expect_false("pltdf" %in% names(r))
  assign("pltdf", NULL, envir = pm_once)
  expect_message(a <- r$pltdf, "pred_df")
  expect_identical(a, r$pred_df)
  expect_silent(r[["pltdf"]])
  p <- ggplot2::autoplot(r)
  expect_match(p$labels$y, "Median of yield")
  expect_match(ggplot2::autoplot(r, scale = "model")$labels$y, "log\\(yield\\)")
  set.seed(2)
  d <- data.frame(A = gl(2, 20), x = c(runif(20, 0, 1), runif(20, 0.5, 2)))
  d$y <- d$x + rnorm(40)
  ra <- q(covariatemeans(lm(y ~ A * x, d), "A", covariate = "x", as_is = TRUE))
  for (l in levels(d$A)) {
    m <- ra$means[ra$means$A == l, ]
    expect_true(all(m$x >= min(d$x[d$A == l]) & m$x <= max(d$x[d$A == l])))
  }
  rf <- q(covariatemeans(lm(y ~ A * x, d), "A", covariate = "x"))
  key <- function(t) paste(t$A, t$x)
  expect_equal(ra$means, rf$means[key(rf$means) %in% key(ra$means), ], ignore_attr = TRUE)
})

test_that("covariatemeans(): two covariates need a matrix 'at'; surface plots", {
  set.seed(1)
  d <- data.frame(A = gl(2, 30), x1 = runif(60), x2 = runif(60))
  d$y <- d$x1 + d$x2 + rnorm(60)
  f <- lm(y ~ A * (x1 + x2), d)
  expect_error(covariatemeans(f, "A", covariate = c("x1", "x2"), at = 1:3), "matrix")
  r <- q(covariatemeans(f, "A", covariate = c("x1", "x2"),
                        at = as.matrix(expand.grid(seq(0, 1, 0.25), seq(0, 1, 0.25)))))
  expect_true(built(ggplot2::autoplot(r)))
  expect_error(ggplot2::autoplot(r, type = "lines"), "one covariate")
})

test_that("semipred(): layouts for two factors, slopes, differences, surfaces, 3d", {
  skip_on_cran()
  set.seed(1)
  dd <- expand.grid(x = seq(0, 1, length = 15), f1 = factor(c("p", "q")), f2 = factor(c("u", "v", "w")), rep = 1:3)
  dd$id <- interaction(dd$f1, dd$f2, dd$rep)
  dd$y <- sin(2*pi*dd$x) + (dd$f1 == "q")*sin(pi*dd$x)^2 + rnorm(nrow(dd), 0, 0.3)
  fx <- q(semireg(y ~ s(x, f1*f2, bs = "os", k = 6) + (1 | id), data = dd))
  r <- q(semipred(fx, "f1:f2", "x"))
  expect_s3_class(ggplot2::autoplot(r)$facet, "FacetWrap")
  expect_true(built(ggplot2::autoplot(r, group = "f2", panel = "f1")))
  expect_true(built(ggplot2::autoplot(r, panel = c("f1", "f2"))))
  s <- q(semipred(fx, "f1", "x", deriv = 1))
  ps <- ggplot2::autoplot(s)
  expect_match(ps$labels$y, "Slope")
  expect_true("GeomHline" %in% vapply(ps$layers, function(l) class(l$geom)[1], ""))
  dfp <- q(semipred(fx, "f1", "x", contr = c("q", "p")))
  expect_match(ggplot2::autoplot(dfp)$labels$y, "Difference")
  expect_error(ggplot2::autoplot(dfp, type = "lines"), "difference")
  set.seed(6)
  d2 <- expand.grid(u = seq(0, 1, length = 7), v = seq(0, 1, length = 7), G = factor(c("a", "b")), r = 1:2)
  d2$w <- sin(3*d2$u) + cos(3*d2$v) + 0.5*(d2$G == "b")*d2$u + rnorm(nrow(d2), sd = 0.2)
  f2 <- q(semireg(w ~ G + s(u, v, by = G, bs = "ztps", k = 5), data = d2))
  sr <- q(semipred(f2, "G", c("u", "v")))
  expect_true(built(ggplot2::autoplot(sr)))
  expect_true(built(ggplot2::autoplot(q(semipred(f2, "G", c("u", "v"), contr = c(2, 1))))))
  skip_if_not_installed("plotly")
  w <- ggplot2::autoplot(sr, type = "3d")
  expect_s3_class(w, "plotly")
  # one scene per level, placed side by side; plotly builds it without warnings
  b <- expect_no_warning(plotly::plotly_build(w))
  expect_setequal(vapply(b$x$data, `[[`, "", "scene"), c("scene", "scene2"))
  expect_equal(b$x$layout$scene2$domain$x, c(0.5, 1))
})

test_that("semipred(): sm_term as block names, labels, terms as written or covariates; advice", {
  skip_on_cran()
  set.seed(2)
  da <- data.frame(x1 = runif(150), x2 = runif(150), z = rnorm(150))
  da$y <- sin(3*da$x1) + da$z*sin(2*pi*da$x2) + rnorm(150, 0, 0.3)
  fa <- q(semireg(y ~ s(x1, bs = "os") + s(x2, by = z, bs = "os"), data = da))
  a <- q(semipred(fa, covariate = "x1", sm_term = "s_x1"))
  for (t in c("s(x1)", "s(x1, bs = \"os\")", "x1")) {
    expect_equal(q(semipred(fa, covariate = "x1", sm_term = t))$means, a$means)
  }
  expect_equal(a$settings$sm_term, "s_x1")
  expect_error(semipred(fa, covariate = "x1", sm_term = "x9"), "s\\(x1\\)")
  expect_false("smooth_other_covariate" %in% advice(a)$id)
  all_sm <- q(semipred(fa, covariate = "x1"))
  expect_true("smooth_other_covariate" %in% advice(all_sm)$id)
  data(Dialyzer, package = "nlme")
  fd <- q(semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (pressure | Subject), data = Dialyzer))
  b <- q(semipred(fd, "QB", "pressure", sm_term = "s(pressure, by = QB)"))
  expect_equal(b$means, q(semipred(fd, "QB", "pressure"))$means)
  expect_equal(q(semipred(fd, "QB", "pressure", sm_term = "s(pressure):QB300"))$settings$sm_term, "s_pressure_QB300")
})

test_that("user layers with other data (only the covariate) can be added to a curve plot", {
  fx <- lm(mpg ~ hp * factor(am), data = mtcars)
  p <- covariatemeans(fx, "factor(am)", covariate = "hp")$plt + ggplot2::geom_rug(data = mtcars, sides = "b")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("marks where the band leaves out the reference are drawn, also on the log axis of ratios", {
  set.seed(7319)
  d <- data.frame(x = rep(seq(0, 1, length = 40), 2), g = factor(rep(c("a", "b"), each = 40)))
  d$cnt <- rpois(80, exp(1 + (d$g == "b") * 2 * d$x))
  r <- covariatemeans(glm(cnt ~ g * x, poisson, d), "g", covariate = "x", contr = c("b", "a"))
  expect_gt(nrow(r$significant), 0)
  p <- ggplot2::autoplot(r)
  k <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomRug"), logical(1)))
  expect_length(k, 1)
  ld <- ggplot2::layer_data(p, k)
  expect_true(nrow(ld) > 0 && all(is.finite(ld$x)))
  # the marks cover the stretches in r$significant
  expect_lte(min(ld$x), min(r$significant$from)); expect_gte(max(ld$x), max(r$significant$to))
})

test_that("autoplot() takes the plot settings of the call; one legend for many groups", {
  skip_on_cran()
  d <- get(utils::data(DataExam5.2, package = "predictmeans", envir = environment()))
  m <- lm(ht ~ sitemean * seedlot, data = d)       # (with site in the model the curves are not estimable)
  cm <- suppressWarnings(covariatemeans(m, "seedlot", "sitemean", plot_args = list(trellis = FALSE, ci = FALSE)))
  p <- ggplot2::autoplot(cm)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], "")
  expect_false("GeomRibbon" %in% geoms)                       # ci = FALSE kept
  expect_s3_class(p$facet, "FacetNull")                       # trellis = FALSE kept
  expect_false(grepl("CI", p$labels$title))                   # no band, no CI in the title
  expect_true("GeomRibbon" %in% vapply(ggplot2::autoplot(cm, interval = "ci")$layers, function(l) class(l$geom)[1], ""))
  # one scale for colour and fill, all 37 seedlots, in the order of the data
  b <- ggplot2::ggplot_build(p)
  sc <- b$plot$scales$get_scales("colour")
  expect_identical(sc$get_limits(), levels(factor(cm$means$seedlot, levels = unique(cm$means$seedlot))))
  expect_length(unique(sc$map(sc$get_limits())), 37)          # no colour repeated
  # the observed points use the same groups (no second set of keys)
  pts <- Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)[[1]]
  expect_identical(levels(pts$data$.group), sc$get_limits())
  # one panel per seedlot: the strips name them, no legend
  p0 <- ggplot2::autoplot(suppressWarnings(covariatemeans(m, "seedlot", "sitemean")))
  expect_identical(p0$theme$legend.position, "none")
})
