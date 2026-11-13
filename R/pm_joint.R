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
#' \item{\code{lmer}}{Kenward-Roger (the default for REML fits): the scaled F test of Kenward and Roger (1997) with the
#'   observed information of the variance parameters, as SAS \code{DDFM=KR}; Satterthwaite (the default for ML fits):
#'   as \code{lmerTest::contest(joint = TRUE)}; containment: the Wald F with the containment df of \code{modelterm}.}
#' \item{KR-expected}{(\code{lmer}, \code{lme}, \code{gls}) Kenward-Roger with the expected information, as
#'   \code{pbkrtest::KRmodcomp()} (see \code{\link{predictmeans}}).}
#' \item{KR-average}{(\code{lmer}, \code{lme}, \code{gls}) Kenward-Roger with the average information: Genstat's REML
#'   df (Genstat's F statistic is the unscaled Wald F).}
#' \item{between-within}{(\code{lmer}, \code{lme}, \code{gls}) the Wald F with the between-within df of
#'   \code{modelterm} (SAS \code{DDFM=BETWITHIN}; see \code{\link{predictmeans}}).}
#' \item{\code{lme}, \code{gls} with \code{df_method = "Satterthwaite"}, \code{"Kenward-Roger"}, \code{"KR-expected"} or \code{"KR-average"}}{as for \code{lmer} (Satterthwaite df of the
#'   F test, lmerTest's rule).}
#' \item{\code{lme}, \code{gls}}{Kenward-Roger for REML fits and Satterthwaite for ML fits (the defaults), as for
#'   \code{lmer}; with \code{df_method = "containment"} the Wald F with \pkg{nlme}'s df of \code{modelterm} (as
#'   \code{anova()}; \code{"containment (nlme)"}) or the residual df (\code{gls}).}
#' \item{other models}{the Wald F with the residual df (normal-based when there is none).}
#' }
#' Contrasts that are linear combinations of others are dropped, so \code{df1} is the rank of the set.
#' Without \code{ctrmatrix}, each set is written as each level against the last (SAS's Type 3 basis): Wald and
#' Kenward-Roger F tests do not depend on the basis, but the Satterthwaite df of a test with several df do (Fai and
#' Cornelius 1996), and with this basis they equal SAS's.
#' \code{emmeans::joint_tests()} gives a slightly different Kenward-Roger test (the unscaled F with its own df).
#' @usage
#' pm_joint(model, modelterm, by = NULL, ctrmatrix = NULL,
#'          df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
#'          weights = c("equal", "proportional", "outer", "cells", "flat"), vcov = NULL)
#' @param model A fitted model, as for \code{\link{predictmeans}}.
#' @param modelterm Name (in "quotes") of the model term whose means are compared, e.g. \code{"nitro:Variety"}.
#' @param by Optional name(s) of factor(s) of \code{modelterm}: test the differences among the means within each of
#'   their levels. Without \code{by} and \code{ctrmatrix}, all differences among the means of \code{modelterm} are
#'   tested.
#' @param ctrmatrix Optional contrast matrix over the means of \code{modelterm} (as for \code{\link{contrastmeans}},
#'   e.g. from \code{\link{contrast_family}} or \code{\link{specify_contrast}}); its rows are tested together.
#' @param df_method For \code{lmer} models (\code{"between-within"} also for \code{lme} and \code{gls}): the
#'   denominator df (see Details).
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
                     df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
                     weights = c("equal", "proportional", "outer", "cells", "flat"), vcov = NULL) {
  df_method <- match.arg(df_method)
  weights <- match.arg(weights)
  if (inherits(model, "aovlist")) model <- pm_aov_carrier(model)
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
    # the same hypotheses, written as each level against the last (SAS's Type 3
    # basis): Wald and Kenward-Roger F do not depend on the basis, but the
    # Satterthwaite df of a several-df test do (Fai and Cornelius 1996)
    sets <- lapply(sets, pm_vs_last)
  }

  ## Model parameters, covariance and df (R/inference.R)
  inf <- pm_inference(model, df_method, vcov_adjust = FALSE, vcov = vcov)
  K <- KK[, inf$mp$estimable, drop = FALSE]
  K[!est_flag, ] <- 0
  # a set with a mean over empty cells: tested on the cells all its means have (R/empty-cells.R)
  over <- rep(NA_character_, length(sets))
  res <- t(vapply(seq_along(sets), function(s) {
    C <- sets[[s]]
    cr <- pm_common_rows(km$empty_cells, C, C %*% K, colnames(K), together = TRUE)
    if (any(cr$sub)) over[s] <<- if (any(cr$na)) "none" else cr$over[which(cr$sub)[1]]
    if (any(cr$na)) return(c(qr(t(C))$rank, NA, NA, NA))
    inf$ftest(cr$L, modelterm)
  }, numeric(4)))
  out <- data.frame(df1 = res[, 1], df2 = res[, 2], F.ratio = res[, 3], p.value = res[, 4])
  if (any(!is.na(over))) {
    out$over <- over          # tested on the cells common to the means (empty cells)
    message("Means over empty cells: tested on the cells all the compared means have (column 'over'; ",
            "\"none\": no cell in common, NA).")
  }
  if (!is.null(by_tab)) out <- cbind(by_tab, out)
  attr(out, "df_method") <- inf$label
  if (!is.null(inf$vcov_label)) attr(out, "vcov") <- inf$vcov_label
  out
}

## Contrasts over means spanning the same space as C, as each level used by C
## against the last of them.
pm_vs_last <- function(C) {
  used <- which(colSums(abs(C)) > 1e-12)
  if (length(used) < 2) return(C)
  k <- length(used)
  M <- matrix(0, k - 1, ncol(C), dimnames = list(colnames(C)[used[-k]], colnames(C)))
  M[cbind(seq_len(k - 1), used[-k])] <- 1
  M[, used[k]] <- -1
  M
}
