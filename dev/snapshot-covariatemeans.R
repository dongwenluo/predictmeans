# Snapshot covariatemeans() output (result tables, plot layer data, plot
# pages) for option combinations written in the 1.1.1 spelling, rewritten to
# the 2.0 spelling when the library has it (as vignette("migrating") says).
# Usage: Rscript --no-init-file dev/snapshot-covariatemeans.R <library> <out.rds>
# then compare two snapshots with all.equal(a, b, tolerance = 0).
lib <- commandArgs(TRUE)[1]; out <- commandArgs(TRUE)[2]
suppressPackageStartupMessages({library(predictmeans, lib.loc = lib); library(ggplot2)})
Oats <- as.data.frame(nlme::Oats)
data(Drug); data(ATP)
m_lm <- lm(yield ~ Block + Variety*nitro, Oats)
m_lmer <- lmerTest::lmer(yield ~ Variety*nitro + (1|Block/Variety), data = Oats)
m_glm <- glm(cbind(R, N - R) ~ Drug + log10Dose, binomial, Drug)
m_log <- lm(log(yield) ~ Variety*nitro, Oats)
Oats$ly1 <- log(Oats$yield + 1)
m_ly1 <- lm(ly1 ~ Variety*nitro, Oats)
m_2x <- lm(yield ~ Block + Variety*nitro, transform(Oats, B2 = factor(as.integer(Block) %% 2)))
calls <- list(
  def = quote(covariatemeans(m_lm, "Variety", covariate = "nitro")),
  noterm = quote(covariatemeans(m_lm, covariate = "nitro")),
  at = quote(covariatemeans(m_lm, "Variety", covariate = "nitro", covariateV = seq(0, 0.6, 0.2))),
  asis = quote(covariatemeans(m_lm, "Variety", covariate = "nitro", as.is = TRUE)),
  level10 = quote(covariatemeans(m_lmer, "Variety", covariate = "nitro", level = 0.1)),
  df20 = quote(covariatemeans(m_lmer, "Variety", covariate = "nitro", Df = 20)),
  glm = quote(covariatemeans(m_glm, "Drug", covariate = "log10Dose")),
  glm_I = quote(covariatemeans(m_glm, "Drug", covariate = "log10Dose", trans = I)),
  log = quote(covariatemeans(m_log, "Variety", covariate = "nitro")),
  off = quote(covariatemeans(m_ly1, "Variety", covariate = "nitro", trans = exp, transOff = 1)),
  plotset = quote(covariatemeans(m_lm, "Variety", covariate = "nitro", trellis = FALSE, mtitle = "T", ci = FALSE,
                                 point = FALSE, plot = TRUE)),
  plotord = quote(covariatemeans(m_lmer, "Variety", covariate = "nitro", plotord = 1)),
  contr = quote(covariatemeans(m_lmer, "Variety", covariate = "nitro", contr = c(1, 2), level = 0.1)),
  contr_sim = quote(covariatemeans(m_glm, "Drug", covariate = "log10Dose", contr = c(1, 2), band = "simultaneous",
                                   contr_scale = "response", mtitle = "D")),
  contr_df = quote(covariatemeans(m_lm, "Variety", covariate = "nitro", contr = c(1, 3), Df = 30)))
pm2 <- function(cl) {
  if (!"conf.level" %in% names(formals(covariatemeans))) return(cl)
  a <- as.list(cl); n <- names(a)
  if ("level" %in% n) { a$conf.level <- 1 - a$level; a$level <- NULL }
  if ("Df" %in% n) { a$df_method <- a$Df; a$Df <- NULL }
  if ("covariateV" %in% n) { a$at <- a$covariateV; a$covariateV <- NULL }
  if ("transOff" %in% n) {
    off <- a$transOff; a$transOff <- NULL
    a$trans <- eval(bquote(function(x) .(a$trans)(x) - .(off)))
  }
  pl <- intersect(n, c("trellis", "plotord", "mtitle", "ci", "point", "jitterv"))
  if (length(pl)) { a$plot_args <- as.call(c(as.name("list"), a[pl])); a[pl] <- NULL }
  as.call(a)
}
plotdata <- function(p) if (inherits(p, "ggplot")) lapply(seq_along(p$layers), function(i) ggplot2::layer_data(p, i)) else p
R <- list()
for (nm in names(calls)) {
  set.seed(42); f <- tempfile(fileext = ".pdf"); grDevices::pdf(f)
  r <- tryCatch(suppressMessages(suppressWarnings(eval(pm2(calls[[nm]])))),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  grDevices::dev.off()
  pages <- length(gregexpr("/Type /Page[^s]", paste(readLines(f, warn = FALSE), collapse = "\n"))[[1]])
  if (is.list(r)) {
    r <- unclass(r)
    r$settings$version <- NULL
    set.seed(44); r$plt <- plotdata(r$plt)
    r$plt_labels <- NULL
  }
  R[[nm]] <- list(result = r, pages = pages)
}
saveRDS(R, out)
cat(length(R), "calls;", sum(vapply(R, function(x) is.character(x$result), logical(1))), "errors\n")
