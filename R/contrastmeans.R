#' Linear Contrast Tests for a Linear Model
#'
#' @description
#' t tests (or permutation t tests) of given contrasts among the predicted means of a model term, for the
#'   models of \code{\link{predictmeans}} (\code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, \code{lmer},
#'   \code{glmer}, \code{glmmTMB}, \code{gam}).
#' @usage
#' contrastmeans(model, modelterm, ctrmatrix, ctrnames=NULL, adj="none", permlist,
#'               vcov_adjust=TRUE,
#'               weights=c("equal", "proportional", "outer", "cells", "flat"),
#'               df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment",
#'       "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
#'               joint=FALSE, conf_level=0.95, vcov=NULL, equivalence=NULL, ...)
#' @param model A fitted model, as for \code{\link{predictmeans}}.
#' @param modelterm Name (in "quotes") for indicating which factor term's contrast to be calculated. 
#'   The \code{modelterm} must be given exactly as it appears in the printed model, e.g. "A" or "A:B".
#' @param ctrmatrix A specified contrast matrix. If it is missing, an edit window asks for it (interactive sessions only). The columns can be named by level (e.g. from \code{\link{specify_contrast}}); levels that are not given get zero weight.
#' @param ctrnames Names of the specified contrasts, e.g. c("A vs D", "C vs B", ...)
#' @param adj Name (in "quote") for indicating a method for adjusting p-values of pairwise comparisons. 
#'   The choices are "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "sidak", "scheffe",
#'   "mvt" and "westfall", in any letter case. "sidak", "scheffe" (rank of \code{ctrmatrix}) and "mvt" (single-step multivariate t, the exact counterpart of SAS
#'   \code{ADJUST=SIMULATE}; needs package \pkg{mvtnorm}) are as in \pkg{emmeans}. With \code{permlist}, "mvt" gives the
#'   single-step and "westfall" the step-down max-|t| adjustment of the permutation p-values (Westfall and Young,
#'   1993); "westfall" needs \code{permlist}, and "scheffe" is not available with it.
#'   The default method is "none".
#' @param permlist A list from \code{\link{permmodels}}: each contrast is then tested by permutation, with p-value \eqn{(b + 1)/(B + 1)} (\eqn{b} of the \eqn{B} permuted |t| at least the observed one), and \code{df_method} is not used.
#'   With the result of \code{permmodels(method = "freedman-lane")} (\code{lm} and
#'   \code{aov} models), each contrast is tested by Freedman-Lane permutations of the model in which that contrast is
#'   0. With a \code{permlist}, \code{vcov_adjust} is not used (the permuted refits use the unadjusted covariance).
#' @param vcov_adjust For models fitted by REML (\code{lmer}, \code{lme}, \code{gls}) with Kenward-Roger df: if \code{TRUE} (the default), standard errors use the Kenward-Roger
#'   adjusted covariance matrix of the fixed effects, as SAS \code{DDFM=KR} and \pkg{emmeans} do, to go with the
#'   Kenward-Roger degrees of freedom. \code{FALSE} uses the unadjusted covariance and reproduces results of
#'   predictmeans 1.x. It has no effect for other models.
#' @param weights Passed to \code{\link{Kmatrix}}: how the cells of factors not in \code{modelterm} are weighted when averaging
#'   over them (\code{"equal"}, the default, \code{"proportional"}, \code{"outer"}, \code{"cells"} or \code{"flat"}).
#' @param df_method Denominator df, by method or as one number (a fixed df, \code{Df} before 2.0): \code{"auto"} (Kenward-Roger for REML fits of \code{lmer}, \code{lme} and \code{gls},
#'   Satterthwaite for ML fits), \code{"Kenward-Roger"}, \code{"Satterthwaite"} or \code{"containment"} (SAS \code{DDFM=CONTAIN}), and
#'   for \code{lmer}, \code{lme} and \code{gls} models \code{"between-within"} (SAS \code{DDFM=BETWITHIN}), and
#'   \code{"KR-expected"} (Kenward-Roger with the expected information, as \pkg{pbkrtest}; \code{"Kenward-Roger"} uses
#'   the observed information, as SAS \code{DDFM=KR}), \code{"KR-average"} (the average information, as Genstat); for \code{lme} and \code{gls} also \code{"Satterthwaite"},
#'   \code{"Kenward-Roger"}, \code{"KR-expected"}, \code{"KR-average"} and \code{"containment"} (\pkg{nlme}'s or the residual df, as
#'   predictmeans 1.x); see \code{\link{predictmeans}}.
#' @param joint Logical: also test all the contrasts together (all of them zero), as \code{\link{pm_joint}} does: an F
#'   test with Kenward-Roger (REML fits, the scaled F of Kenward and Roger, 1997) or Satterthwaite (ML fits)
#'   denominator df, the residual df for other models (normal-based, i.e. chi-square /
#'   df1, when there is none). Contrasts that are linear combinations of others are dropped, so \code{df1} is the rank
#'   of \code{ctrmatrix}. With \code{permlist}, also a permutation p-value of the joint Wald F statistic, from the same
#'   permutations (for Freedman-Lane, the model with all the contrasts zero). The default is \code{FALSE}.
#' @param conf_level Confidence level of the limits of the contrasts (\code{level = 0.05} before 2.0); the tests
#'   are at the significance level \code{1 - conf_level}. The limits are t limits of each contrast, or with \code{permlist} permutation-t limits:
#'   estimate +- c SE, with c the critical value of the permutation test (the k-th largest permuted |t|,
#'   k = floor((1 - conf_level) (nperm + 1))). With \code{adj} "bonferroni", "sidak", "scheffe" or "mvt" (with \code{permlist}:
#'   "bonferroni", "sidak", "mvt" or "westfall", from the largest permuted |t| of the contrasts) the table also has
#'   \code{adj.lower.CL} and \code{adj.upper.CL}: simultaneous limits that exclude 0 exactly when the adjusted
#'   p-value is below \code{1 - conf_level}, as \pkg{emmeans}' \code{confint()} and SAS \code{AdjLower}/\code{AdjUpper}. The
#'   step-down and false discovery rate methods have none.
#' @param vcov Optional covariance of the coefficients to use instead of the model's. For \code{lm}, \code{aov}
#'   and \code{glm} models: \code{"HC0"}, \code{"HC1"}, \code{"HC2"} or \code{"HC3"} (heteroscedasticity-consistent,
#'   as \code{sandwich::vcovHC()}), or a one-sided formula naming a cluster variable of the data, e.g. \code{~ Block}
#'   (cluster-robust CR1, as \code{sandwich::vcovCL()}). For any model: a covariance matrix of the coefficients, or
#'   a function that returns one from the model (for mixed models e.g. from the \pkg{clubSandwich} package). The
#'   degrees of freedom stay those of the model (the residual df; for \code{lmer} Kenward-Roger or Satterthwaite df
#'   computed from the model), and \code{vcov_adjust} is not used. With a \code{permlist}, \code{vcov} is not used
#'   (for robust permutation tests see \code{permmodels(robust = TRUE)}). The default \code{NULL} uses the model's
#'   covariance.
#' @param equivalence Optional equivalence margins for equivalence tests (two one-sided t tests, TOST) of the
#'   contrasts: a positive number \eqn{d} (equivalence region \eqn{-d} to \eqn{d}) or \code{c(lower, upper)}, on the
#'   scale of the comparisons (the link scale for GLMs: e.g. log ratios). The TOST p-value is the larger of the two
#'   one-sided p-values; below \code{1 - conf_level}, the difference lies within the margins (as emmeans'
#'   \code{test(..., side = "equivalence")}). Not adjusted for multiplicity; with a \code{permlist} it is t-based
#'   (normal-based in \code{contrastmeans}).
#' @param ... Not used: arguments removed in 2.0 (\code{Df}, \code{level}) stop with their 2.0 spelling.
#' @return
#' An object of class \code{"contrastmeans"}, a list. \code{print()} shows the contrasts, \code{tidy()} and
#' \code{as.data.frame()} return the table \code{contrasts} (\code{what = "joint"}: the joint test). Its first
#' components are those of 1.x:
#' \item{The t tests of the specified contrasts}{A table (rounded to 4 decimals) of the t tests of the contrasts,
#'   with the adjusted p-value in \code{Pr(>|t|)}.}
#' \item{K}{A contrast matrix.}
#' \item{Joint F test of the contrasts}{With \code{joint = TRUE}: \code{df1}, \code{df2}, \code{F.ratio} and
#'   \code{p.value} of the joint test (and \code{Perm_p_value} with \code{permlist}).}
#' and then the tables in the layout of \code{\link{predictmeans}}, unrounded:
#' \item{contrasts}{One row per contrast: \code{contrast} (its name), \code{over} (contrasts of means over empty
#'   cells: the cells used), \code{estimate}, \code{std.error}, \code{df}, \code{statistic} (t),
#'   \code{p.value} (unadjusted; with \code{permlist} the permutation p-value), \code{p.adj} (adjusted by
#'   \code{adj}), \code{p.equiv} (with \code{equivalence}), \code{conf.low}, \code{conf.high},
#'   \code{adj.conf.low}, \code{adj.conf.high} (simultaneous, with the adjustments that have them),
#'   \code{p.method}, and for a log or logit scale (link, or a log response) \code{ratio} or
#'   \code{odds.ratio} with its limits when the weights of every contrast sum to zero.}
#' \item{joint}{With \code{joint = TRUE}: the joint test, unrounded.}
#' \item{settings}{\code{modelterm}, \code{level} (significance level, \code{1 - conf_level}), \code{adjust}, \code{df_method}, \code{vcov_adjust}, \code{vcov},
#'   \code{equivalence}, \code{nperm}, \code{model_class}, \code{version}.}
#' @references
#' Torsten Hothorn, Frank Bretz and Peter Westfall (2008),
#'     \emph{Simultaneous Inference in General Parametric Models. Biometrical}, Journal 50(3), 346--363.
#' @seealso \code{\link{specify_contrast}}, \code{\link{contrast_family}}, \code{\link{pm_joint}}, \code{\link{predictmeans}}, \code{\link{tidy.predictmeans}}, \code{\link{advice}}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#'
#' # 1. Enter a contrast matrix in an edit window (interactive sessions only)
#' if (interactive()) contrastmeans(fm, "nitro:Variety")
#'  
#' # 2. Give the contrast matrix: Victory - Golden Rain at nitro 0, Victory at 0 - Marvellous at 0.4
#' cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#'             c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#' contrastmeans(fm, "nitro:Variety", ctrmatrix=cm)
#'
#' # 3. ... and both contrasts tested together (an F test with 2 numerator df)
#' ct <- contrastmeans(fm, "nitro:Variety", ctrmatrix=cm, joint=TRUE)
#' ct
#' as.data.frame(ct)
#' @name contrastmeans
NULL

#' @export
contrastmeans <- function(model, modelterm, ctrmatrix, ctrnames=NULL, adj="none", permlist, vcov_adjust=TRUE,
                          weights=c("equal", "proportional", "outer", "cells", "flat"),
                          df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"), joint=FALSE,
                          conf_level=0.95, vcov=NULL, equivalence=NULL, ...) {
  pm_defunct("contrastmeans", match.call(), strict = TRUE, formals_of = formals(sys.function()))
  level <- pm_conf_alpha(conf_level)
  vcov_adjust_default <- missing(vcov_adjust)
  dfa <- pm_df_arg(df_method, eval(formals(sys.function())$df_method))
  df_method <- dfa$method
  df_user <- !is.null(dfa$Df)
  if (df_user) Df <- dfa$Df
  pm_equiv_margins(equivalence)                  # check early
  pm_note_singular(model)
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)

  model <- pm_means_model(model)

  if (is.null(ctrnames) || all(ctrnames%in%c("NULL", ""))) {
    ctrnames <- NULL
  }
  weights <- match.arg(weights)
  adj <- pm_match_adj(adj, tukey = FALSE, westfall = TRUE)
  if (adj == "westfall" && missing(permlist)) {
    stop("adj = \"westfall\" (step-down max-|t|) needs permutation p-values ('permlist'); without them use ",
         "\"mvt\" (single step) or \"holm\".", call. = FALSE)
  }
  km <- Kmatrix(model, modelterm, weights=weights)
  KK <- km$K
  termsLabel <- rownames(KK)
  if (missing(ctrmatrix)) {
    if (!interactive()) {
      stop("Please supply 'ctrmatrix'; the contrast editor needs an interactive session.")
    }
    ques <- paste("\nHow many contrasts you want set up for ", sQuote(modelterm), "? ", sep="")
    nctr <- as.integer(readline(ques))
    nrK <- length(termsLabel)
    xnew <- matrix(rep(0, nctr*nrK), ncol=nrK)
    colnames(xnew) <- termsLabel
    xnew <- edit(xnew)
    ctrmatrix <- xnew[1:nctr, 1:nrK, drop=FALSE]
    rownames(ctrmatrix) <- ctrnames
    if (!(all(rowSums(ctrmatrix)==0)))  {
      cat("\n", "The contrast matrix is:\n\n")
      print(ctrmatrix)
      cat("\n\n")
      stop("\n", "Please check the row",  sQuote(which(rowSums(ctrmatrix)!=0)), "of the contrast matrix!\n\n")
    }
  } else {
    if (is.null(dim(ctrmatrix))) {
      ctrmatrix <- matrix(ctrmatrix, nrow=1, dimnames=list(NULL, names(ctrmatrix)))
    }
    if (!is.null(colnames(ctrmatrix)) && all(colnames(ctrmatrix) %in% termsLabel)) {
      # columns named by level (e.g. from specify_contrast()): put them in the
      # order of the means, with zero weight for levels not mentioned
      full <- matrix(0, nrow(ctrmatrix), length(termsLabel), dimnames=list(rownames(ctrmatrix), termsLabel))
      full[, colnames(ctrmatrix)] <- ctrmatrix
      ctrmatrix <- full
    } else {
      if (ncol(ctrmatrix) != length(termsLabel)) {
        stop("'ctrmatrix' must have one column for each of the ", length(termsLabel), " levels of ",
             sQuote(modelterm), ": ", paste(termsLabel, collapse=", "))
      }
      colnames(ctrmatrix) <- termsLabel
    }
    if (!is.null(ctrnames)) rownames(ctrmatrix) <- ctrnames
    nctr <- nrow(ctrmatrix)
  }
  if (all((is.null(dim(ctrmatrix)) | dim(ctrmatrix)[1]==1), missing(permlist)))  {
    adj <- "none"
  }

  # Means for missing cells (NaN rows of K) or other non-estimable means can
  # only be left out of a contrast, not used in it.
  est_flag <- attr(KK, "estimable")
  if (is.null(est_flag)) est_flag <- rep(TRUE, nrow(KK))
  bad_row <- !est_flag | apply(KK, 1, function(r) any(is.na(r)))
  used_bad <- colSums(abs(ctrmatrix[, bad_row, drop=FALSE])) > 0
  if (any(used_bad)) {
    stop("The contrast(s) use mean(s) that cannot be estimated: ",
         paste(termsLabel[bad_row][used_bad], collapse=", "))
  }
  KK[bad_row, ] <- 0
  inf <- pm_inference(model, df_method, vcov_adjust, vcov_adjust_default, vcov,
                      permlist = if (!missing(permlist)) permlist)          # R/inference.R
  mp <- inf$mp
  K <- pm_K_cols(KK, mp)
  rK <- ctrmatrix%*%K
  # contrasts with a mean over empty cells: on the cells all their means have (R/empty-cells.R)
  cr <- pm_common_rows(km$empty_cells, ctrmatrix, rK, colnames(K))
  rK <- cr$L
  if (any(cr$sub)) {
    cn <- if (is.null(rownames(ctrmatrix))) paste("contrast", seq_len(nrow(ctrmatrix))) else rownames(ctrmatrix)
    message("Contrasts with means over empty cells use only the cells their means all have: ",
            paste0(cn[cr$sub], ": ", ifelse(cr$na[cr$sub], "none (NA)", cr$over[cr$sub]), collapse = "; "))
  }

  if (missing(permlist) && !df_user) Df <- inf$df(rK, modelterm)

  cm <- rK%*%mp$coef
  vcovm <- mp$vcov
  vcov.contr <- rK %*% tcrossprod(vcovm, rK)
  ses <- sqrt(diag(vcov.contr))
  cm[cr$na] <- NA                 # no cells in common
  ses[cr$na] <- NA
  t.v <- cm/ses
  nr <- nrow(rK)
  dv <- t(1/ses)
  cor.contr <- as.matrix(vcov.contr * (t(dv) %*% dv))

  if (missing(permlist)) {
   # t.p.value <- 2*pt(-abs(t.v), Df)
	t.p.value <- apply(cbind(t.v, Df), 1, function(x) 2*pt(-abs(x[1]), x[2]))
    p_raw <- t.p.value
    # one family (R/p-adjust.R); Scheffe with the rank of ctrmatrix
    fams <- pm_families(rep("all", nr), est = !is.na(as.vector(cm)), rank = qr(ctrmatrix)$rank)
    t.p.value <- pm_adjust_p(adj, t.v, t.p.value, rep_len(Df, nr), fams, R = cor.contr)
    crit <- stats::qt(1 - level/2, Df)
    # simultaneous limits that agree with the adjusted p-values (R/p-adjust.R)
    crit_adj <- pm_adjust_crit(adj, level, rep_len(Df, nr), fams, R = cor.contr)
    eqm <- pm_equiv_margins(equivalence)
    p_equiv <- if (!is.null(eqm)) pm_tost(as.vector(cm), ses, Df, eqm)
    # the table, once; the 1.x table is read from it
    tab <- pm_contrast_table(rownames(ctrmatrix), cm, ses, rep_len(Df, nr), t.v, p_raw, t.p.value, crit, ses, crit_adj,
                             p_equiv, "t", cr, pm_contrast_kind(model, ctrmatrix))
    out.put <- pm_contrast_legacy(tab, rownames(ctrmatrix), perm = FALSE)
    attr(out.put,"Note") <- paste("The p-value is adjusted by", sQuote(adj), "method, if p-value = 0 means p-value < 0.0001.")
    if (!is.null(inf$vcov_label)) {
      attr(out.put, "Note") <- paste0(attr(out.put, "Note"), " Standard errors from the covariance given by 'vcov' (",
                                      inf$vcov_label, ").")
    }
    }else{
  	nsim <- pm_perm_n(permlist)
  if (adj == "scheffe") stop("adj = \"scheffe\" is not available for permutation p-values; use another 'adj'.",
                             call. = FALSE)
  # permutation p-values; for adj = "mvt" / "westfall" the joint permutation
  # distribution of the contrasts' |t| (max-|t|, Westfall-Young)
  maxt <- nrow(rK) > 1 && adj %in% pm_perm_maxt_methods
  fam <- pm_perm_family(permlist, rK, model, abs(as.vector(t.v)), 1, joint = TRUE)
  per.p <- fam$p
  p_raw <- per.p
  fams <- pm_families(rep("all", nr), est = !is.na(as.vector(cm)))
  if (nrow(rK) > 1 && adj != "none") {
    per.p <- pm_adjust_p(adj, t.v, per.p, rep(NA_real_, nr), fams, perm_T = list(obs = fam$obs, T = fam$T))
  }
  # permutation-t limits of each contrast, and simultaneous ones with the
  # adjustments that have them (max-|t|, bonferroni, sidak; R/p-adjust.R)
  crit <- pm_perm_crit(fam$T, level)
  ci_se <- if (!is.null(fam$se)) fam$se else ses   # robust permlist: HC3 SEs
  crit_adj <- if (nrow(rK) > 1) pm_adjust_crit(adj, level, rep(NA_real_, nr), fams, T = fam$T)
  eqm <- pm_equiv_margins(equivalence)
  p_equiv <- if (!is.null(eqm)) pm_tost(as.vector(cm), ses, Inf, eqm)       # normal-based
  tab <- pm_contrast_table(rownames(ctrmatrix), cm, ses, rep(NA_real_, nr), t.v, p_raw, per.p, crit, ci_se, crit_adj,
                           p_equiv, if (identical(attr(permlist, "resampling"), "parametric bootstrap")) "bootstrap" else
                             "permutation", cr, pm_contrast_kind(model, ctrmatrix))
  out.put <- pm_contrast_legacy(tab, rownames(ctrmatrix), perm = TRUE)
  attr(out.put,"Note") <- paste0("The permuted p-value is obtained using ", sQuote(nsim), " permutations",
                                 if (nrow(rK) > 1 && adj != "none") paste0(", adjusted by ", sQuote(adj), " method"),
                                 if (maxt) paste0(" (", if (adj == "westfall") "step-down" else "single-step",
                                                  " max-|t| of the permutations)"), ".")
  }
  out <- list("The t tests of the specified contrasts"=round(out.put, 4), "K"=ctrmatrix)
  jt_raw <- NULL
  if (isTRUE(joint)) {
    # all the contrasts together: the F test of pm_joint(), and with a
    # permlist the permutation p-value of the joint Wald F
    jt <- pm_joint(model, modelterm, ctrmatrix = ctrmatrix, df_method = df_method, weights = weights, vcov = vcov)
    if (!missing(permlist)) jt$Perm_p_value <- pm_perm_p(fam$Fobs, fam$F, fam$exact)
    dfm <- attr(jt, "df_method")
    jt_raw <- jt
    rownames(jt_raw) <- NULL
    jt[] <- lapply(jt, function(v) if (is.numeric(v)) round(v, 4) else v)
    rownames(jt) <- "all contrasts"
    attr(jt, "Note") <- paste0("F test of all the contrasts together (", nrow(ctrmatrix), " contrasts, rank ", jt$df1,
                               "; denominator df: ", dfm,
                               if (!is.null(attr(jt, "vcov"))) paste0("; covariance: ", attr(jt, "vcov")), ")",
                               if (!missing(permlist)) paste0("; Perm_p_value from ", sQuote(nsim), " permutations"),
                               ".")
    out[["Joint F test of the contrasts"]] <- jt
  }
  # the stable tables (the 1.x elements above keep their places)
  out$contrasts <- tab
  out$joint <- jt_raw
  out$settings <- list(modelterm = modelterm, level = level, adjust = adj,
                       df_method = if (!missing(permlist)) "permutation" else if (df_user) "user" else inf$label,
                       vcov_adjust = isTRUE(inf$vcov_adjust),
                       vcov = if (is.null(inf$vcov_label)) "model" else inf$vcov_label,
                       equivalence = eqm, nperm = if (!missing(permlist)) nsim else NA_integer_,
                       model_class = class(model)[1], version = as.character(utils::packageVersion("predictmeans")))
  class(out) <- "contrastmeans"
  out <- pm_add_advice(out, pm_contrast_advice(model, tab, adj))
  out
}

## The table of contrastmeans() in the columns of predictmeans()$pairs:
## contrast, over (contrasts on common cells), estimate, std.error, df,
## statistic, p.value (unadjusted), p.adj, p.equiv, conf.low, conf.high
## (each contrast's limits), adj.conf.low, adj.conf.high (simultaneous, with
## the adjustments that have them), p.method, and for a log or logit scale
## the ratio or odds ratio with its limits (contrasts whose weights sum to
## zero).
pm_contrast_table <- function(names, est, se, df, t, p, p_adj, crit, ci_se, crit_adj, p_equiv, p_method, cr, kind) {
  est <- as.vector(est)
  n <- length(est)
  out <- data.frame(contrast = if (is.null(names)) paste("contrast", seq_len(n)) else names, stringsAsFactors = FALSE)
  if (any(cr$sub)) out$over <- ifelse(cr$sub, cr$over, NA_character_)
  out$estimate <- est
  out$std.error <- as.vector(se)
  out$df <- as.vector(df)
  out$statistic <- as.vector(t)
  out$p.value <- as.vector(p)
  out$p.adj <- as.vector(p_adj)
  if (!is.null(p_equiv)) out$p.equiv <- p_equiv
  out$conf.low <- est - crit*ci_se
  out$conf.high <- est + crit*ci_se
  if (!is.null(crit_adj)) {
    out$adj.conf.low <- est - crit_adj*ci_se
    out$adj.conf.high <- est + crit_adj*ci_se
  }
  out$p.method <- p_method
  if (!is.null(kind)) {
    pre <- if (kind == "log") "ratio" else "odds.ratio"
    out[[pre]] <- exp(est)
    out[[paste0(pre, ".conf.low")]] <- exp(out$conf.low)
    out[[paste0(pre, ".conf.high")]] <- exp(out$conf.high)
  }
  out
}

## The 1.x table of contrastmeans() (before rounding), read from the table of
## pm_contrast_table(): the adjusted p-value in "Pr(>|t|)" ("Permuted
## Pr(>|t|)" with a permlist, which has no df column).
pm_contrast_legacy <- function(tab, rn, perm = FALSE) {
  cols <- c(Estimate = "estimate", "Std. Error" = "std.error", "t value" = "statistic",
            if (!perm) c(df = "df"), stats::setNames("p.adj", if (perm) "Permuted Pr(>|t|)" else "Pr(>|t|)"),
            lower.CL = "conf.low", upper.CL = "conf.high", adj.lower.CL = "adj.conf.low",
            adj.upper.CL = "adj.conf.high", p.equiv = "p.equiv")
  cols <- cols[cols %in% names(tab)]
  m <- as.matrix(tab[cols])
  dimnames(m) <- list(rn, names(cols))
  m
}

## "log" or "logit" when exp of the contrasts is a ratio or an odds ratio: a
## log or logit link, or a log response in the formula, and every contrast's
## weights summing to zero. NULL otherwise.
pm_contrast_kind <- function(model, ctrmatrix) {
  if (any(abs(rowSums(ctrmatrix)) > 1e-8)) return(NULL)
  link_scale <- !is.null(pm_link_trans(model))
  trans <- if (link_scale) NULL else pm_response_trans(model)$fun
  pm_ratio_kind(trans, model, link_scale)
}

## Methods: the contrasts (or the joint test) as a data frame, and a printed
## summary.
#' @exportS3Method generics::tidy contrastmeans
tidy.contrastmeans <- function(x, what = c("contrasts", "joint"), ...) {
  pm_result_table(x, if (!missing(what)) what, c("contrasts", "joint"))    # R/tidy.R
}

#' @exportS3Method as.data.frame contrastmeans
as.data.frame.contrastmeans <- function(x, row.names = NULL, optional = FALSE, what = c("contrasts", "joint"), ...) {
  pm_result_table(x, if (!missing(what)) what, c("contrasts", "joint"), row.names)
}

#' @exportS3Method print contrastmeans
print.contrastmeans <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("Contrasts of '", s$modelterm, "' (", format(100*(1 - s$level)), "% limits; p.adj: ", s$adjust, "; ",
      if (identical(s$df_method, "permutation")) paste(s$nperm, "permutations") else paste(s$df_method, "df"), ")\n",
      sep = "")
  tab <- x$contrasts
  tab$p.method <- NULL
  print(tab, digits = digits, row.names = FALSE)
  jt <- x$joint
  if (!is.null(jt)) {
    cat("Joint F test: F(", jt$df1, ", ", format(round(jt$df2, 1)), ") = ", format(round(jt$F.ratio, 3)), ", p = ",
        format.pval(jt$p.value, digits = 3), if (!is.null(jt$Perm_p_value)) paste0("; permutation p = ",
                                                                                     format.pval(jt$Perm_p_value, digits = 3)),
        "\n", sep = "")
  }
  if (!identical(s$vcov, "model")) cat("Standard errors from the covariance given by 'vcov' (", s$vcov, ").\n", sep = "")
  pm_print_advice(x)
  invisible(x)
}
