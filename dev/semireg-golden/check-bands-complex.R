## Contrast bands of semipred() for GLMMs with correlated random effects:
## semipred(contr = ) against the covariance from lme4's own pieces in the
## spherical form (b = Lambda u, u ~ N(0, I)): V(beta, u) = (Cu'W Cu +
## diag(0, I))^-1 with Cu = [X, Z Lambda], mapped to (beta, b); this needs no
## inverse of Lambda Lambda' (singular when a correlation is +-1). Poisson
## counts, three treatments' curves over time, 8 blocks with their own
## intercepts and slopes, a crossed day effect.
## Result (2026-11-17, unconditional = FALSE):
##   (a) (t | Block)               difference exact, SE 6.5e-16 (correlation estimated at 1)
##   (b) (s(t) | Block) curves     difference exact, SE 3.2e-16
##   (c) (t | Block) + (1 | Day)   difference exact, SE 2.3e-16
##   (d) semireg_tmb() on (a): logLik equal, difference 8e-4, SE 1e-4 (optimiser tolerance)
##   Rscript dev/semireg-golden/check-bands-complex.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
q <- function(e) suppressWarnings(suppressMessages(e))
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
set.seed(31)
dd <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:8), Obs = factor(1:2))
dd$Day <- factor(dd$t %% 4)                                  # crossed with Block
bi <- rnorm(8, 0, 0.3); bs <- rnorm(8, 0, 0.04)
eta <- 1 + sin(dd$t / 2) / 2 + 0.04 * (dd$A == "a2") * dd$t + bi[dd$Block] + bs[dd$Block] * dd$t +
  rnorm(4, 0, 0.2)[dd$Day]
dd$y <- rpois(nrow(dd), exp(eta))
tg <- seq(0, 10, by = 0.5)
## semipred(contr) against V = (C'WC + (Lambda Lambda')^-1)^-1 from lme4's pieces
check <- function(fp, label) {
  m <- fp$semer
  pr <- q(semipred(fp, "A", "t", covariateV = tg, contr = c("a2", "a1"), unconditional = FALSE, prt = FALSE))$pred_df
  Xm <- lme4::getME(m, "X"); Zm <- lme4::getME(m, "Z"); b <- lme4::getME(m, "b")
  etaf <- as.vector(Xm %*% lme4::fixef(m) + Zm %*% b)
  w <- exp(etaf)                                             # Poisson, log link: mu.eta^2 / V = mu
  kx <- ncol(Xm)
  Lam <- as.matrix(lme4::getME(m, "Lambda"))
  Cu <- cbind(Xm, as.matrix(Zm) %*% Lam)                  # (beta, u), b = Lambda u, u ~ N(0, I)
  Vu <- solve(crossprod(Cu * w, Cu) + diag(c(rep(0, kx), rep(1, ncol(Lam)))))
  Tm <- as.matrix(Matrix::bdiag(diag(kx), Lam))
  V <- Tm %*% Vu %*% t(Tm)
  C <- cbind(Xm, as.matrix(Zm))
  LL <- tcrossprod(Lam)
  des <- q(pm_semi_design(fp, "A", "t", tg, NULL, fp$data, unconditional = FALSE))
  D <- as.matrix(des$D_mx); g <- des$pred_df$A
  gp <- lme4::getME(m, "Gp"); cn <- names(lme4::getME(m, "cnms"))
  cols <- unlist(lapply(des$sm_term, function(bk) { i <- match(bk, cn); kx + (gp[i] + 1):gp[i + 1] }))
  L <- matrix(0, sum(g == "a2"), ncol(C)); L[, c(seq_len(kx), cols)] <- D[g == "a2", ] - D[g == "a1", ]
  est <- as.vector(L %*% c(lme4::fixef(m), as.vector(b))); se <- sqrt(rowSums((L %*% V) * L))
  offd <- sum(abs(LL[upper.tri(LL)]) > 0)
  cat(sprintf("%-38s off-diagonal entries of Lambda Lambda': %4d | difference %.1e, SE %.1e | df %s\n", label, offd,
              rel(pr$diff, est), rel(pr$ses, se), paste(unique(round((pr$UL - pr$diff) / pr$ses, 3)), collapse = ",")))
  vc <- as.data.frame(lme4::VarCorr(m)); vc <- vc[!grepl("^s_", vc$grp) | grepl("Block", vc$grp), ]; cat("     ", paste(sprintf("%s %s%s %.3g", vc$grp, vc$var1, ifelse(is.na(vc$var2), "", paste0(",", vc$var2)), vc$sdcor), collapse = "; "), "\n")
  invisible(pr)
}
fa <- q(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (t | Block), family = poisson, data = dd, prt = FALSE))
pa <- check(fa, "(a) (t | Block)")
fb <- q(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (s(t, bs = "os", k = 4) | Block), family = poisson, data = dd, prt = FALSE))
check(fb, "(b) (s(t) | Block)")
fc <- q(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (t | Block) + (1 | Day), family = poisson, data = dd, prt = FALSE))
check(fc, "(c) (t | Block) + (1 | Day)")
## (d) the glmmTMB engine on (a)
ft <- q(semireg_tmb(y ~ s(t, by = A, bs = "os", k = 6) + (t | Block), family = poisson, data = dd, prt = FALSE))
pt <- q(semipred(ft, "A", "t", covariateV = tg, contr = c("a2", "a1"), unconditional = FALSE, prt = FALSE))$pred_df
cat(sprintf("(d) glmmTMB engine vs lme4 engine, (t | Block): logLik %.4f vs %.4f | difference %.1e, SE %.1e\n",
            as.numeric(logLik(ft$semer)), as.numeric(logLik(fa$semer)), rel(pt$diff, pa$diff), rel(pt$ses, pa$ses)))
