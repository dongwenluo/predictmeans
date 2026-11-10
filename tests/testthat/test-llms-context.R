# The context file for AI assistants must list every column of the stable tables.
test_that("inst/llms-context.md documents every column of means and pairs", {
  f <- system.file("llms-context.md", package = "predictmeans")
  skip_if(!nzchar(f), "llms-context.md not found")
  doc <- paste(readLines(f, warn = FALSE), collapse = "\n")
  documented <- function(cols) cols[!vapply(cols, function(x) grepl(paste0("`", x, "`"), doc, fixed = TRUE), logical(1))]
  data(Drug, envir = environment())
  g <- glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial, data = Drug)
  pg <- predictmeans(g, "Drug", pairwise = TRUE, plot = FALSE, prtnum = FALSE)
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  fm <- lm(log(yield) ~ Block + nitro * Variety, data = d)
  pl <- predictmeans(fm, "nitro:Variety", atvar = "Variety", trans = exp, plot = FALSE, prtnum = FALSE)
  factor_cols <- c("Drug", "log10Dose", "nitro", "Variety")
  expect_length(documented(setdiff(union(names(pg$means), names(pl$means)), factor_cols)), 0)
  expect_length(documented(union(names(pg$pairs), names(pl$pairs))), 0)
  expect_length(documented(names(pg$settings)), 0)
})
