## The table of arguments and functions removed in 2.0 (R/defunct.R).

test_that("the defunct table is well formed", {
  tab <- predictmeans:::pm_defunct_table
  expect_false(anyDuplicated(tab[c("fun", "arg")]) > 0)
  ns <- asNamespace("predictmeans")
  expect_true(all(vapply(unique(tab$fun), exists, logical(1), envir = ns, inherits = FALSE)))
  # a value is filled in only where the 2.0 spelling has a place for it
  expect_true(all(grepl("%s", tab$use[nzchar(tab$value)], fixed = TRUE)))
  expect_true(all(tab$value %in% c("", "same", "one_minus")))
})

test_that("a removed argument stops with its 2.0 spelling and the value given", {
  d <- predictmeans:::pm_defunct
  expect_error(d("predictmeans", quote(predictmeans(m, "A", level = 0.05))),
               "'level' was removed in predictmeans 2.0: use conf_level = 0.95", fixed = TRUE)
  expect_error(d("predictmeans", quote(predictmeans(m, "A", level = 0.1))), "conf_level = 0.9", fixed = TRUE)
  expect_error(d("contrastmeans", quote(contrastmeans(m, "A", Df = 10))), "use df_method = 10", fixed = TRUE)
  k <- 12
  expect_error(d("covariatemeans", quote(covariatemeans(m, covariateV = list(x = k))), env = environment()),
               "use at = list(x = k)", fixed = TRUE)
  expect_error(d("permlmer", quote(permlmer(lmer0 = a, lmer1 = b))), "model0 = a.*\n.*model1 = b")
  expect_error(d("semireg", quote(semireg(y ~ x, smoothZ = list()))), "smooths in the formula")
  expect_error(d("predictmeans", quote(predictmeans(m, "A", newwd = TRUE))), "vignette(\"migrating\")", fixed = TRUE)
  # kept arguments and other functions pass
  expect_null(d("predictmeans", quote(predictmeans(m, "A", adj = "tukey"))))
  expect_null(d("pm_check", quote(pm_check(m, level = 0.05))))
})

test_that("strict: an unknown name stops as an unused argument", {
  f <- function(model, plot = FALSE, ...) predictmeans:::pm_defunct("pm_check", match.call(), strict = TRUE,
                                                                     formals_of = formals(sys.function()))
  expect_null(f(1, plot = TRUE))
  expect_error(f(1, plt = TRUE), "unused argument: plt")
})

test_that("a removed function stops with its replacement", {
  for (fn in names(predictmeans:::pm_defunct_funs)) {
    expect_error(predictmeans:::pm_defunct_fun(fn), paste0(fn, "() was removed"), fixed = TRUE, class = "defunctError")
  }
})

test_that("no removed argument can be partially matched to a remaining one", {
  # Once an argument is no longer a formal, the function must take '...' (so
  # that pm_defunct() sees it), and it must not be the start of a formal
  # placed before '...' (R would pass it there silently).
  tab <- predictmeans:::pm_defunct_table
  ns <- asNamespace("predictmeans")
  for (i in seq_len(nrow(tab))) {
    fm <- names(formals(get(tab$fun[i], envir = ns)))
    if (tab$arg[i] %in% fm) next
    expect_true("..." %in% fm, info = paste(tab$fun[i], tab$arg[i]))
    before <- fm[seq_len(match("...", fm) - 1L)]
    expect_false(any(startsWith(before, tab$arg[i])), info = paste(tab$fun[i], tab$arg[i]))
  }
})
