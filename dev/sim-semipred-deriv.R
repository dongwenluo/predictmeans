## Coverage of the bands of derivatives of fitted curves, semipred(deriv = ):
## one gaussian curve, y = sin(x/1.5) + N(0, 0.3^2), 150 points with x
## uniform on [0, 10] (both ends included), semireg(y ~ s(x, bs = "os", k)).
## True slope cos(x/1.5)/1.5, curvature -sin(x/1.5)/2.25; 37 points from 0.5
## to 9.5. Pointwise coverage (mean over the points) and simultaneous
## coverage (the band covers at every point), for the O'Sullivan basis with
## 6 interior knots (the default) and with 15.
##   Rscript dev/sim-semipred-deriv.R [nsim] [out.rds]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
nsim <- if (length(args)) as.integer(args[1]) else 500
q <- function(e) suppressWarnings(suppressMessages(e))
xs <- seq(0.5, 9.5, by = 0.25)
truth <- list(function(x) cos(x/1.5)/1.5, function(x) -sin(x/1.5)/2.25)
one <- function(seed, k) {
  set.seed(seed)
  d <- data.frame(x = c(0, 10, runif(148, 0, 10)))
  d$y <- sin(d$x/1.5) + rnorm(150, 0, 0.3)
  f <- q(semireg(y ~ s(x, bs = "os", k = k), data = d, prt = FALSE))
  out <- NULL
  for (dv in 1:2) for (band in c("pointwise", "simultaneous")) {
    p <- q(semipred(f, covariate = "x", covariateV = xs, deriv = dv, band = band, prt = FALSE))$pred_df
    inside <- p$LL <= truth[[dv]](xs) & p$UL >= truth[[dv]](xs)
    out <- c(out, if (band == "pointwise") mean(inside) else all(inside))
  }
  stats::setNames(out, c("d1_pw", "d1_sim", "d2_pw", "d2_sim"))
}
res <- list()
for (k in c(6, 15)) {
  r <- parallel::mclapply(seq_len(nsim), function(i) tryCatch(one(4000 + i, k), error = function(e) NULL),
                          mc.cores = parallel::detectCores())
  r <- do.call(rbind, r[!vapply(r, is.null, logical(1))])
  res[[paste0("k", k)]] <- r
  m <- colMeans(r)
  cat(sprintf("k = %2d (%d data sets): slope pointwise %.3f simultaneous %.3f | curvature pointwise %.3f simultaneous %.3f\n",
              k, nrow(r), m[["d1_pw"]], m[["d1_sim"]], m[["d2_pw"]], m[["d2_sim"]]))
}
if (length(args) > 1) saveRDS(res, args[2])
