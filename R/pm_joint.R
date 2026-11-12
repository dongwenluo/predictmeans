#' Joint F Tests of Contrasts among Means (Slices)
#'
#' @description
#' F tests of sets of contrasts among the predicted means of a model term. With \code{by}, the means of the term are
#' compared within each level of the \code{by} factor(s): the test of \code{nitro} within each \code{Variety}, as SAS
#' \code{LSMEANS nitro*Variety / SLICE=Variety} or \code{emmeans::joint_tests(model, by = "Variety")}. With
#' \code{ctrmatrix}, all its rows are tested together.
#' @details
#' The F statistic and its denominator degrees of freedom follow the model:
#' \describe{
#' \item{\code{lmer}}{Kenward-Roger (the default for REML fits): the scaled F test of Kenward and Roger (1997), as
#'   \code{pbkrtest::KRmodcomp()} and the ANOVA tables of this package; Satterthwaite (the default for ML fits): as
#'   \code{lmerTest::contest(joint = TRUE)}; containment: the Wald F with the containment df of \code{modelterm}.}
#' \item{\code{lme}}{the Wald F with the smallest \pkg{nlme} df of \code{modelterm} and its factors (as
#'   \code{\link{contrastmeans}}).}
#' \item{other models}{the Wald F with the residual df (normal-based when there is none).}
#' }
#' Contrasts that are linear combinations of others are dropped, so \code{df1} is the rank of the set.
#' \code{emmeans::joint_tests()} gives a slightly different Kenward-Roger test (the unscaled F with its own df).
#' @usage
#' pm_joint(model, modelterm, by = NULL, ctrmatrix = NULL,
#'          df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment"),
#'          weights = c("equal", "proportional", "outer", "cells", "flat"), vcov = NULL)
#' @param model A fitted model, as for \code{\link{predictmeans}}.
#' @param modelterm Name (in "quotes") of the model term whose means are compared, e.g. \code{"nitro:Variety"}.
#' @param by Optional name(s) of factor(s) of \code{modelterm}: test the differences among the means within each of
#'   their levels. Without \code{by} and \code{ctrmatrix}, all differences among the means of \code{modelterm} are
#'   tested.
#' @param ctrmatrix Optional contrast matrix over the means of \code{modelterm} (as for \code{\link{contrastmeans}},
#'   e.g. from \code{\link{contrast_family}} or \code{\link{specify_contrast}}); its rows are tested together.
#' @param df_method For \code{lmer} models: the denominator df (see Details).
#' @param weights How the means average over other factors, as in \code{\link{predictmeans}}.
#' @param vcov Optional covariance of the coefficients, as in \code{\link{predictmeans}}: the F statistic is then
#'   the Wald F with this covariance, and the denominator df those of the model's own test.
#' @return
#' A data frame with one row per test: the \code{by} levels (if any), \code{df1}, \code{df2}, \code{F.ratio} and
#' \code{p.value}. Attribute \code{"df_method"} says how \code{df2} was found.
#' @references
#' Kenward, M. G. and Roger, J. H. (1997), Small sample inference for fixed effects from restricted maximum likelihood,
#' \emph{Biometrics}, 53, 983-997.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   # nitrogen within each variety (SAS SLICE)
#'   pm_joint(fm, "nitro:Variety", by = "Variety")
#'   # varieties within each nitrogen level
#'   pm_joint(fm, "nitro:Variety", by = "nitro")
#'   # a set of contrasts: linear and quadratic nitrogen trends together
#'   pm_joint(fm, "nitro", ctrmatrix = contrast_family(Kmatrix(fm, "nitro"), "poly", max.degree = 2))
#' @name pm_joint
NULL

pm_joint <- function(model, modelterm, by = NULL, ctrmatrix = NULL,
                     df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment"),
                     weights = c("equal", "proportional", "outer", "cells", "flat"), vcov = NULL) {
  df_method <- match.arg(df_method)
  weights <- match.arg(weights)
  if (inherits(model, "aovlist")) model <- aovlist_lmer(model)
  pm_check_modelterm(model, modelterm)
  km <- Kmatrix(model, modelterm, weights = weights)
  KK <- km$K
  levs <- rownames(KK)
  est_flag <- attr(KK, "estimable")
  if (is.null(est_flag)) est_flag <- rep(TRUE, nrow(KK))
  est_flag <- est_flag & !apply(KK, 1, function(r) any(is.na(r)))

  ## The sets of contrasts (over the means), one per test
  if (!is.null(ctrmatrix)) {
    if (!is.null(by)) stop("Give either 'by' or 'ctrmatrix', not both.")
    if (is.null(dim(ctrmatrix))) ctrmatrix <- matrix(ctrmatrix, nrow = 1, dimnames = list(NULL, names(ctrmatrix)))
    if (!is.null(colnames(ctrmatrix)) && all(colnames(ctrmatrix) %in% levs)) {
      full <- matrix(0, nrow(ctrmatrix), length(levs), dimnames = list(rownames(ctrmatrix), levs))
      full[, colnames(ctrmatrix)] <- ctrmatrix
      ctrmatrix <- full
    } else if (ncol(ctrmatrix) != length(levs)) {
      stop("'ctrmatrix' must have one column for each of the ", length(levs), " levels of ", sQuote(modelterm), ": ",
           paste(levs, collapse = ", "))
    }
    if (any(colSums(abs(ctrmatrix[, !est_flag, drop = FALSE])) > 0)) {
      stop("The contrast(s) use mean(s) that cannot be estimated: ", paste(levs[!est_flag], collapse = ", "))
    }
    sets <- list(ctrmatrix)
    by_tab <- NULL
  } else {
    fam <- suppressMessages(contrast_family(km, "consec", by = by))
    if (is.null(by)) {
      sets <- list(fam)
      by_tab <- NULL
    } else {
      fac <- km$fctnames
      gkey <- do.call(paste, c(lapply(fac[by], as.character), sep = ":"))
      gname <- unique(gkey)
      grp_of_row <- sub(": .*$", "", rownames(fam))
      gname <- gname[gname %in% grp_of_row]
      sets <- lapply(gname, function(g) fam[grp_of_row == g, , drop = FALSE])
      by_tab <- unique(fac[match(gname, gkey), by, drop = FALSE])
      rownames(by_tab) <- NULL
    }
  }

  ## Model parameters and the df machinery
  dfm <- pm_df_method(model, df_method)
  df_type <- dfm$df_type
  mp <- mymodelparm(model, df = 0)
  Vuser <- pm_vcov_resolve(model, vcov, names(mp$coef))
  K <- KK[, mp$estimable, drop = FALSE]
  K[!est_flag, ] <- 0
  other_df <- NULL
  if (inherits(model, "lme")) {
    vars <- c(unlist(strsplit(modelterm, "\\:")), modelterm)
    other_df <- min(stats::terms(model$fixDF)[vars], na.rm = TRUE)
  } else if (!inherits(model, "lmerMod")) {
    rd <- mymodelparm(model)$df
    other_df <- if (is.numeric(rd) && length(rd) == 1 && is.finite(rd) && rd > 0) rd else Inf
  }
  cont_df <- if (identical(df_type, "containment")) pm_containment_df(model)[[modelterm]] else NULL

  res <- t(vapply(sets, function(C) {
    pm_joint_ftest(model, C %*% K, mp, df_type, dfm$vcov_adj, other_df, cont_df, Vuser)
  }, numeric(4)))
  out <- data.frame(df1 = res[, 1], df2 = res[, 2], F.ratio = res[, 3], p.value = res[, 4])
  if (!is.null(by_tab)) out <- cbind(by_tab, out)
  attr(out, "df_method") <- if (inherits(model, "lmerMod")) df_type else if (inherits(model, "lme")) "nlme" else "residual"
  if (!is.null(Vuser)) attr(out, "vcov") <- attr(Vuser, "label")
  out
}

## F test of L beta = 0: c(df1, df2, F, p). Rows of L that are linear
## combinations of others are dropped first.
pm_joint_ftest <- function(model, L, mp, df_type, vcov_adj, other_df, cont_df, Vuser = NULL) {
  beta <- mp$coef
  V <- as.matrix(mp$vcov)
  r <- qr(t(L))
  q <- r$rank
  if (q == 0) return(c(0, NA, NA, NA))
  L <- L[r$pivot[seq_len(q)], , drop = FALSE]
  if (!is.null(Vuser)) {
    # a given covariance ('vcov'): its Wald F, with the denominator df of the
    # model's own test
    out <- pm_joint_ftest(model, L, mp, df_type, vcov_adj, other_df, cont_df)
    b <- L %*% beta
    Fv <- as.numeric(t(b) %*% solve(L %*% Vuser %*% t(L), b)) / q
    return(c(q, out[2], Fv, stats::pf(Fv, q, out[2], lower.tail = FALSE)))
  }
  wald <- function(Vb) {
    b <- L %*% beta
    as.numeric(t(b) %*% solve(L %*% Vb %*% t(L), b))
  }
  if (inherits(model, "lmerMod") && identical(df_type, "Kenward-Roger")) {
    kr <- pm_kr_ftest(model, L, PhiA = vcov_adj)
    return(unname(c(q, kr["ddf"], kr["stat"], kr["p.value"])))
  }
  if (inherits(model, "lmerMod") && identical(df_type, "Satterthwaite")) {
    # lmerTest's rule for several contrasts (contestMD): split L V L' into
    # its eigen-directions, each with a Satterthwaite df nu_m, and combine
    # E = sum(nu_m / (nu_m - 2)), df2 = 2E / (E - q)
    parts <- vcov_adj
    VL <- L %*% parts$vcov_beta %*% t(L)
    eg <- eigen(VL, symmetric = TRUE)
    PtL <- crossprod(eg$vectors, L)
    Fv <- sum(drop(PtL %*% beta)^2/eg$values)/q
    nu <- pm_satt_ddf(PtL, parts)
    E <- sum((nu/(nu - 2))[nu > 2])
    ddf <- if (E > q) 2*E/(E - q) else NA_real_
    return(c(q, ddf, Fv, stats::pf(Fv, q, ddf, lower.tail = FALSE)))
  }
  Fv <- wald(V)/q
  ddf <- if (!is.null(cont_df)) cont_df else other_df
  if (is.null(ddf)) ddf <- Inf
  c(q, ddf, Fv, stats::pf(Fv, q, ddf, lower.tail = FALSE))
}
