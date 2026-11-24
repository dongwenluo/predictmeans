test_that("every id in the code is listed on the help page ?advice", {
  src <- list.files(test_path("..", "..", "R"), full.names = TRUE)
  skip_if(!length(src), "package sources not available")
  code <- unlist(lapply(src, readLines, warn = FALSE))
  ids <- unique(regmatches(code, regexpr('(?<=pm_advice\\(")[a-z0-9_]+', code, perl = TRUE)))
  rd <- paste(readLines(test_path("..", "..", "man", "advice.Rd"), warn = FALSE), collapse = "\n")
  missing <- ids[!vapply(ids, function(i) grepl(paste0("code{", i, "}"), rd, fixed = TRUE), logical(1))]
  expect_identical(missing, character(0))
  expect_identical(advice(1), pm_advice_get(1))
})

test_that("pm_advice() makes one row and checks its input", {
  a <- pm_advice("x", "Some text.")
  expect_equal(a, data.frame(id = "x", level = "note", text = "Some text."))
  expect_equal(pm_advice("y", "t", "caution")$level, "caution")
  expect_error(pm_advice("y", "t", "other"))
  expect_error(pm_advice("", "t"))
  expect_error(pm_advice("y", character(0)))
})

test_that("pm_add_advice() stores, replaces by id and orders cautions first", {
  x <- pm_add_advice(1, pm_advice("a", "note a"), NULL, pm_advice("b", "caution b", "caution"))
  a <- pm_advice_get(x)
  expect_equal(a$id, c("b", "a"))
  x <- pm_add_advice(x, pm_advice("a", "new a"))
  expect_equal(pm_advice_get(x)$text[pm_advice_get(x)$id == "a"], "new a")
  expect_equal(nrow(pm_advice_get(x)), 2)
  expect_identical(pm_add_advice(2), 2)
  expect_equal(nrow(pm_advice_get(3)), 0)
})

test_that("pm_print_advice() prints one block and respects the option", {
  x <- pm_add_advice(1, pm_advice("b", "Caution text.", "caution"), pm_advice("a", "Note text."))
  out <- capture.output(pm_print_advice(x))
  expect_equal(out[2], "Advice:")
  expect_equal(out[3], "! Caution text.")
  expect_equal(out[4], "i Note text.")
  op <- options(predictmeans.advice = FALSE)
  res <- tryCatch(capture.output(pm_print_advice(x)), finally = options(op))
  expect_length(res, 0)
  expect_length(capture.output(pm_print_advice(1)), 0)
})
