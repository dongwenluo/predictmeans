#' Matrix of Coefficients in a Linear Model
#'
#' @description
#' This function obtains a matrix of coefficients (the linear functions \eqn{K} such that the
#'   predicted means are \eqn{K\hat\beta}) for models such as \code{aov}, \code{lm},
#'   \code{glm}, \code{gls}, \code{lme}, \code{lmer}, \code{glmer}, \code{glmmTMB}, \code{gam}
#'   (package \pkg{mgcv}) and \code{semireg}.
#'
#'   The predicted means are built on a reference grid: every combination of the factor levels in
#'   the model, with covariates fixed at their means (or other chosen values). Each row of \code{K}
#'   is the average of the grid rows (weighted as \code{weights} says; equal by default) for one level combination of
#'   \code{modelterm}. Combinations of nested factors that cannot exist (a term of the model holds both
#'   factors, e.g. \code{G/T}, and each level of one occurs with one level of the other) are left out of
#'   the averages; so are, with \code{weights = "equal"}, the empty cells of interaction terms (see
#'   \code{\link{predictmeans}}). When a mean cannot be estimated because two factors that share no term of the
#'   model are confounded in the data (not all their combinations occur: a factor nested in another only in the
#'   data, each pasture on one farm), the combinations of those two factors without data are left out as well,
#'   with a warning that names them: the mean is the average over the combinations that exist. The pair is found
#'   by trying, pairs among the factors averaged over first. Writing the nesting in the formula (\code{G/T}) or
#'   leaving out the aliased factor gives the same means; \code{weights = "cells"} weights the combinations by
#'   their numbers of observations. The design matrix of the grid is evaluated with the model's own stored basis
#'   information, so covariates may enter the model through transformations such as
#'   \code{poly(x, 2)}, \code{ns(x, 3)}, \code{bs(x)}, \code{log(x)}, \code{scale(x)} or
#'   \code{I(x^2)}, only through interactions (e.g. \code{A:x} without \code{x}, or \code{x:z}),
#'   or through smooth terms (\code{s(x)}, \code{te(x, z)}, \code{s(x, by = A)}) of a \code{gam}.
#' @usage
#' Kmatrix(model, modelterm, covariate=NULL, at=NULL, data=NULL, cov_reduce=mean,
#'         exclude=NULL, weights=c("equal", "proportional", "outer", "cells", "flat"),
#'         refgrid=FALSE, ...)
#' @param model A fitted model: \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, \code{lmer},
#' \code{glmer}, \code{glmmTMB}, \code{gam}/\code{bam} (package \pkg{mgcv}; for \code{gamm} use its
#' \code{$gam} component) or \code{semireg}.
#' @param modelterm Name (in "quotes") of the variables whose combinations give the rows of \code{K}: fixed-effect
#'   variables of the model joined by ":", e.g. "A" or "A:B". It need not be a term of the model: any combination
#'   of its factors can be asked for (e.g. "A:B" for \code{y ~ A + B}), and a covariate may be among them (its
#'   values as given in \code{covariate}, else its mean).
#'   A factor created in the model formula, e.g. \code{factor(dose)}, can be given either as
#'   "factor(dose)" or as "dose".
#' @param covariate The values of covariates; those not mentioned are at their means (see \code{cov_reduce}).
#'   \itemize{
#'   \item A named list (or named numeric vector) of values, e.g. \code{list(x = c(-1, 0, 1), z = 0.3)}: a covariate
#'     with one value is fixed there; one with several values gives rows of \code{K} for each value (all
#'     combinations for several such covariates). A data frame gives joint points: one set of rows for each of
#'     its rows.
#'   \item An unnamed numeric vector with one value for each covariate of the model, in the order they have in
#'     the model. With another length this is an error (before 2.0 the values were recycled, and names were
#'     ignored).
#'   \item The name of one or more covariates, in which case predictions are made over a range of values of those
#'     covariates (see \code{at}).
#'   }
#' For a \code{semireg} model this may also be a covariate that only appears in the smooth terms.
#' @param at Values for the covariate(s) named in \code{covariate} (\code{covariateV} before 2.0). For one covariate, a numeric vector. For several covariates,
#'   either a list of numeric vectors (all combinations are used; named by covariate, or in the order of
#'   \code{covariate}) or a matrix with one named column per covariate,
#'   where each row is one prediction point (e.g. a set of locations for a two-dimensional smooth).
#'   When \code{at} is \code{NULL}, 50 (one covariate) or 30 (several covariates) equally spaced values over the observed range are used.
#' @param data In some cases, you need to provide the data set used in model fitting, for example when the data
#' used to fit the model can no longer be found. Usually the data are recovered from the model.
#' @param ... Not used: arguments removed in 2.0 (\code{covariateV}, \code{prtnum}) stop with their 2.0 spelling.
#' @param cov_reduce How covariates that are not varied are handled. A function returning a single value, applied to each covariate
#' (the default \code{mean}; \code{median} is another choice), or the string \code{"average"}, which averages the
#' predictions over the observed values of those covariates instead of fixing them at one value. With
#' \code{"average"}, nonlinear terms such as \code{poly(x, 2)} or \code{s(x)} are averaged over the data, which for a
#' \code{gam} gives predicted means adjusted for (rather than evaluated at the centre of) a spatial or temporal trend.
#' @param exclude For \code{gam} models only: labels of smooth terms to leave out of the predictions (their columns of \code{K} are
#' set to zero), e.g. \code{"s(Block)"}. The default \code{NULL} excludes random-effect smooths (\code{bs = "re"} and
#' \code{bs = "fs"}), so the predicted means are population-level means, as for \code{lmer}. Use
#' \code{exclude = character(0)} to keep all smooth terms.
#' @param weights How the cells of the factors that are not in \code{modelterm} are weighted when averaging over them:
#'   \code{"equal"} (the default; equal weights), \code{"proportional"} (weights proportional to the frequencies of
#'   those cells in the data, the same for every mean, as SAS \code{LSMEANS / OM}), \code{"outer"} (products of each
#'   factor's own frequencies), \code{"cells"} (each mean weights the cells by their own frequencies) or
#'   \code{"flat"} (equal weights over the non-empty cells). The definitions are those of \pkg{emmeans}.
#' @param refgrid Return the grid of all combinations (element \code{refgrid})? The default \code{FALSE} builds
#'   it only when it is needed (\code{semireg} models, means over empty cells). Otherwise each group of variables
#'   that share a model term is averaged on its own grid, which gives the same \code{K} and is much smaller
#'   when factors do not interact: for \code{y ~ Line + Block + Site} with 300, 90 and 4 levels, grids of 300,
#'   90 and 4 rows instead of 108,000.
#' @return
#' \item{K}{Coefficients matrix. It has an attribute \code{"estimable"}, a logical vector with one value per row that is
#'  \code{FALSE} for rows that are not estimable functions of the model parameters (e.g. means for a missing cell in a
#'  rank-deficient model).}
#'  \item{fctnames}{A model frame contains factor(s) info in the model.}
#'  \item{response}{The name of response variable in the model.}
#'  \item{preddf}{The data of the predictor variables used in model fitting.}
#'  \item{refgrid}{\code{NULL}, or (see argument \code{refgrid}) the reference grid of all combinations: a list with the grid data (\code{grid}), the sparse
#'  weight matrix \code{W} that averages grid rows into rows of \code{K}, and the rows of \code{fctnames}
#'  (\code{target}), and \code{basis}, the function that gives the model's basis for grid rows. Any other basis
#'  evaluated on \code{grid}, e.g. the smooth terms of a \code{semireg} model, can be averaged in the same way as
#'  \code{W \%*\% basis}.}
#'  \item{empty_cells}{\code{NULL}, or when some means average over empty cells (\code{weights = "equal"}; see
#'  \code{\link{predictmeans}}): \code{affected}, which rows of \code{K} do so, \code{missing}, the empty cells, and
#'  \code{compare(C, together = FALSE)}, which gives comparisons \code{C} (rows over the rows of \code{K}) on the
#'  cells their means have in common (\code{together = TRUE}: the cells common to all the means of \code{C}).}
#' @references
#' Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004),
#' 	\emph{Prediction in linear mixed models}, Australian and New Zealand Journal of Statistics, 46(3), 325-347.
#' @seealso \code{\link{predictmeans}}, \code{\link{contrast_family}}, \code{\link{specify_contrast}}
#' @author
#' This function heavily depends on the codes from package "lsmeans".
#' @examples
#'   library(predictmeans)
#'   data(Oats, package="nlme")
#'   fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#'   Kmatrix(fm, "Variety")$K
#'   Kmatrix(fm, "Variety", 0.5)$K
#'   Kmatrix(fm, "Variety", "nitro")$K
#'   Kmatrix(fm, "Variety", "nitro", at=seq(0, 0.6, 0.1))$K
#'   Kmatrix(fm, "Variety", list(nitro = c(0, 0.3, 0.6)))$K      # the same, values by name
#'   # not a term of the model
#'   Kmatrix(lm(yield ~ Block + Variety + nitro, data=Oats), "Block:Variety")$fctnames
#'
#'   # Covariates entering through transformations or only through interactions
#'   fm2 <- lm(yield ~ Block + Variety + poly(nitro, 2), data=Oats)
#'   Kmatrix(fm2, "Variety")$K
#'   fm3 <- lm(yield ~ Variety + Variety:nitro, data=Oats)
#'   Kmatrix(fm3, "Variety", "nitro", at=c(0, 0.3, 0.6))$K
#'
#'   # A smooth term in an mgcv model, averaged over the observed nitro values
#'   if (requireNamespace("mgcv", quietly = TRUE)) {
#'     fm4 <- mgcv::gam(yield ~ Variety + s(nitro, k = 4) + s(Block, bs = "re"), data = Oats)
#'     Kmatrix(fm4, "Variety", cov_reduce = "average")$K
#'   }
#' @name Kmatrix
NULL

#' @export
Kmatrix <- function(model, modelterm, covariate=NULL, at=NULL, data=NULL, cov_reduce=mean, exclude=NULL,
                    weights=c("equal", "proportional", "outer", "cells", "flat"), refgrid=FALSE, ...) {
  pm_defunct("Kmatrix", match.call(), strict = TRUE, formals_of = formals(sys.function()))

  weights <- match.arg(weights)
  model <- pm_means_model(model)

  ## 1. The model and the data of its fixed-effect variables
  info <- pm_model_info(model, data, exclude)
  model <- info$model
  Terms <- info$Terms
  rec <- pm_recover_data(info)
  model_data <- rec$data
  pv <- pm_predictors(info, model_data, rec$needed)
  X <- pv$X

  ## 2. The rows of K: the variables of 'modelterm' (factors: their levels;
  ## covariates: their values) and the values of the other covariates
  term_labels <- unlist(strsplit(modelterm, "\\:"))
  term_vars <- pm_term_vars(term_labels)
  pm_check_term_vars(term_labels, term_vars, pv)
  average_cov <- identical(cov_reduce, "average")
  if (!average_cov && !is.function(cov_reduce)) {
    stop("'cov_reduce' must be a function such as mean or median, or \"average\".")
  }
  cv <- pm_cov_values(covariate, at, pv$covariates, X)
  baselevs <- pm_base_values(pv, cv, cov_reduce)
  facs <- c(term_vars, setdiff(cv$varied, term_vars))

  # Optionally average over the observed values of the covariates that are
  # not varied or given, instead of fixing them at one value.
  profiles <- prof_w <- NULL
  avg_vars <- character(0)
  if (average_cov) {
    avg_vars <- setdiff(pv$covariates, c(facs, names(cv$fixed)))
    if (length(avg_vars)) {
      pdat <- X[avg_vars]
      pkey <- do.call(paste, c(lapply(pdat, function(z) format(z, digits = 15)), sep = "\r"))
      first <- !duplicated(pkey)
      profiles <- pdat[first, , drop = FALSE]
      rownames(profiles) <- NULL
      prof_w <- as.vector(table(factor(pkey, levels = pkey[first]))) / length(pkey)
    }
  }
  ## 3. The weights of the cells averaged over
  factor_names <- pv$factors
  by_factor <- any(term_vars %in% factor_names)
  # Nested factors: combinations that cannot exist are left out of the
  # averages. The formula says which factors are nested (a term with both, as
  # G/T or G:T). Nesting that is only in the data (y ~ droplet + vol, each vol
  # with one droplet) is confounding: see the fallback after assemble().
  nested_data <- pm_nested_pairs(X, factor_names)
  declared <- vapply(nested_data, function(p) any(vapply(pm_term_sets(Terms), function(s) all(c(p$a, p$b) %in% s),
                                                         logical(1))), logical(1))
  nested <- nested_data[declared]
  cw <- pm_cell_weights(weights, X, factor_names, setdiff(factor_names, facs))
  Terms_rule <- Terms

  ## 4. The basis of the fixed effects for rows of predictor values,
  ## through pm_basis() (see ?pm_coef), with the covariates' values beside it
  basis <- pm_grid_basis(info, model_data, pv$xlev)
  Terms <- basis$Terms
  basis <- basis$fun
  cov_report <- setdiff(pv$covariates, facs)
  lc <- NULL
  basis_cov <- function(nd) {
    if (length(lc)) nd <- pm_set_level_covariates(nd, lc)
    B <- basis(nd)
    if (length(cov_report)) B <- cbind(B, matrix(unlist(lapply(nd[cov_report], as.numeric)), nrow(nd),
                                                 dimnames = list(NULL, paste0("\rcov:", cov_report))))
    B
  }
  dims <- pm_grid_dims(baselevs, cv$joint, profiles, prof_w)
  # the null space of the design matrix of the fit (rank-deficient fits only:
  # aliased or dropped coefficients), for the estimability check
  null_space <- NULL
  if (info$builder != "gam" && info$rank_deficient(colnames(basis(model_data[1, , drop = FALSE])))) {
    null_space <- pm_null_space(basis(model_data))
  }
  # Covariates that are a function of a factor (a site-level value such as
  # the site mean): in a rank-deficient fit they are aliased with that factor,
  # and a grid row that pairs a level with any other value of the covariate
  # is not estimable. Each row then takes its level's value (pm_level_covariates()).
  lc <- if (!is.null(null_space) && !average_cov) {
    pm_level_covariates(X, setdiff(pv$covariates, c(facs, names(cv$fixed))), factor_names)
  }
  if (length(lc)) {
    message("Covariate", if (length(lc) > 1) "s", " ", paste(names(lc), collapse = ", "), " ",
            if (length(lc) > 1) "are" else "is", " constant within the levels of ",
            paste(unique(vapply(lc, `[[`, "", "factor")), collapse = ", "),
            ": each mean uses the value of each level averaged over.")
    basis0 <- basis
    basis <- function(nd) basis0(pm_set_level_covariates(nd, lc))
  }

  ## 5. K: the basis averaged over the reference grid, and 6. estimability.
  ## 'confounded': pairs of factors whose combinations without data are left
  ## out of the averages as empty cells are (see below).
  assemble <- function(confounded = list()) {
    # empty cells of the interaction terms (R/empty-cells.R): with equal
    # weights, means average over the observed cells only
    empty_rule <- if (weights == "equal") pm_empty_rule(Terms_rule, factor_names, X, facs, nested, apply_rule = by_factor,
                                                        confounded = confounded)
    present <- NULL
    if (length(nested) || !is.null(empty_rule) || !is.null(cw$fun)) {
      present <- function(grid) {
        w <- if (is.null(cw$fun)) rep(1, nrow(grid)) else as.numeric(cw$fun(grid))
        if (length(nested)) w <- w * !pm_absent_by_design(grid, nested)
        if (!is.null(empty_rule)) w <- w * empty_rule$present(grid)
        w
      }
    }
    # The grid of all combinations is built only when it is needed
    # ('refgrid', semireg fits, comparisons on common cells); otherwise each
    # group of interacting variables is averaged on its own grid
    # (pm_average_by_groups()).
    need_grid <- isTRUE(refgrid) || isTRUE(info$semireg) || length(empty_rule$missing) > 0 || length(lc) > 0
    groups <- if (!need_grid) pm_dim_groups(dims, c(pm_basis_sets(info, Terms), cw$sets,
                                                   lapply(nested, function(p) c(p$a, p$b))))
    rg <- NULL
    if (need_grid || length(groups) < 2) {
      rg <- pm_make_refgrid(facs, dims, pv$matdat, present)
      Kc <- as.matrix(pm_refgrid_average(rg, basis_cov))
      target <- rg$target
      empty <- rg$empty
    } else {
      av <- pm_average_by_groups(facs, dims, pv$matdat, present, groups, pm_ref_positions(dims, X), basis_cov)
      Kc <- av$K
      target <- av$target
      empty <- av$empty
    }
    is_cov <- startsWith(colnames(Kc), "\rcov:")
    K <- Kc[, !is_cov, drop = FALSE]
    # a row of K is estimable if it is orthogonal to the null space
    estimable <- rep(TRUE, nrow(K))
    null_proj <- NULL
    if (!is.null(null_space)) {
      null_proj <- K[, rownames(null_space), drop = FALSE] %*% null_space
      attr(null_proj, "scale") <- attr(null_space, "scale")
      estimable <- pm_estimable_rows(null_proj, K)
    }
    estimable[is.na(estimable) | empty] <- FALSE
    list(empty_rule = empty_rule, rg = rg, Kc = Kc, K = K, is_cov = is_cov, target = target, empty = empty,
         estimable = estimable, null_proj = null_proj, n_na = sum(!estimable & !empty))
  }
  out <- assemble()
  # Means that cannot be estimated because factors that do not interact in
  # the model are confounded in the data (not all their combinations occur:
  # each pasture on one farm; a spray volume used with one droplet size
  # only): the mean over all their combinations does not exist. As for empty
  # cells, the mean is then the average over the combinations that exist,
  # with a warning that names the factors (the maintainer's decision of
  # 2026-11-16). The pairs are found by trying, in the order below, and the
  # first choice that leaves no mean without an estimate is used.
  if (out$n_na > 0 && weights == "equal" && by_factor) {
    cand <- pm_confounded_pairs(X, factor_names, Terms_rule, nested)
    # first the pairs among the factors averaged over (the mean keeps its
    # meaning: only the averaging changes), then those with a factor of the
    # term (its levels then stand for different levels of the other factor)
    own <- vapply(cand, function(p) any(p %in% facs), logical(1))
    nu <- cand[!own]
    ow <- cand[own]
    tries <- c(lapply(nu, list), if (length(nu) > 1) list(nu), lapply(ow, list),
               if (length(ow) && length(cand) > 1) list(c(nu, ow)))
    best <- out
    for (tr in tries) {
      o <- tryCatch(assemble(tr), error = function(e) NULL)
      if (!is.null(o) && o$n_na < best$n_na) {
        best <- o
        if (o$n_na == 0) break
      }
    }
    out <- best
  }
  empty_rule <- out$empty_rule; rg <- out$rg; Kc <- out$Kc; K <- out$K; is_cov <- out$is_cov
  target <- out$target; empty <- out$empty; estimable <- out$estimable; null_proj <- out$null_proj

  fctnames <- target
  names(fctnames)[seq_along(term_vars)] <- term_labels
  rnK <- do.call("paste", c(fctnames[, seq_along(facs), drop = FALSE], sep=":"))
  rownames(K) <- rnK
  miss <- setdiff(info$coef_names, colnames(K))
  if (length(miss) != 0) {
    stop("The model's coefficients ", paste(utils::head(miss, 5), collapse = ", "), if (length(miss) > 5) ", ...",
         " have no column in the design of '", modelterm, "'. If a factor has levels without data, drop them ",
         "(droplevels()) and refit; for a new model class, pm_coef() and pm_basis() must give the same column ",
         "names (see ?pm_coef).", call. = FALSE)
  }
  names(estimable) <- rnK
  if (any(!estimable & !empty)) pm_note_nesting(nested_data[!declared], facs, weights)

  # means that average over empty cells: say which, and prepare their
  # comparisons on common cells (R/empty-cells.R)
  empty_cells <- if (length(empty_rule$missing)) pm_empty_cells(rg, empty_rule, basis, info$coef_names, null_space, rnK,
                                                                modelterm)

  K <- K[, info$coef_names, drop = FALSE]
  attr(K, "estimable") <- estimable
  # components along non-estimable directions: a contrast of the rows is
  # estimable only when these cancel (used for pairwise differences)
  if (!is.null(null_proj)) {
    rownames(null_proj) <- rnK
    attr(K, "null_proj") <- null_proj
  }

  ## 7. Covariates that are fixed (not varied) are reported with the value used.
  if (length(cov_report)) {
    cov_vals <- Kc[, is_cov, drop = FALSE]
    if (any(empty)) cov_vals[empty, ] <- NaN
    colnames(cov_vals) <- cov_report
    rownames(cov_vals) <- NULL
    fctnames <- cbind(fctnames, cov_vals)
  }
  rownames(fctnames) <- NULL

  if (!is.null(rg)) rg$basis <- basis    # the basis of the grid rows (for means averaged on the response scale)
  return(list(K=K, fctnames=fctnames, response=info$yname, preddf=X, refgrid=rg, data=model_data,
              empty_cells=empty_cells, level_covariates=lc))
}

## ---- The steps of Kmatrix() ----

## The fixed-effect variables of the model, from the data used in the fit:
##   X           one column per variable (characters and logicals as factors)
##   xlev        levels of the factors
##   coerced     numeric variables used as factors in the formula (factor(dose))
##   factors     the variables whose levels are averaged over or define rows
##   covariates  the numeric variables (not offsets)
##   matdat      matrix covariates: their column means
##   excluded_only  variables only in excluded smooth terms (gam)
pm_predictors <- function(info, model_data, needed) {
  form <- as.formula(paste("~", paste(pm_bq(needed), collapse = "+")), env = info$envir)
  X <- model.frame(form, model_data, na.action = na.omit, drop.unused.levels = TRUE)
  for (v in names(X)) {
    if (is.character(X[[v]]) || is.logical(X[[v]])) X[[v]] <- factor(X[[v]])
  }
  coerced <- intersect(pm_coerced_vars(formula(info$Terms)), names(X))
  excluded_only <- intersect(info$excluded_only, names(X))
  is_fac <- vapply(X, is.factor, logical(1))
  is_mat <- vapply(X, is.matrix, logical(1))
  xlev <- lapply(X[is_fac], levels)
  matdat <- lapply(X[is_mat], function(obj) apply(obj, 2, mean, na.rm = TRUE))
  list(X = X, xlev = xlev, coerced = coerced, matdat = matdat, excluded_only = excluded_only,
       factors = setdiff(c(names(xlev), coerced), excluded_only),
       covariates = setdiff(names(X)[!is_fac & !is_mat], c(coerced, info$offset_vars, excluded_only)))
}

## The variables of 'modelterm' must be fixed-effect variables of the model.
pm_check_term_vars <- function(term_labels, term_vars, pv) {
  known <- c(pv$factors, pv$covariates)
  bad <- !term_vars %in% known
  if (any(term_vars %in% pv$excluded_only)) {
    stop(paste("Factor(s)", paste(term_labels[term_vars %in% pv$excluded_only], collapse = ", "),
               "only appear in smooth terms excluded from prediction (see 'exclude')."))
  }
  if (any(bad)) {
    hit <- known[match(tolower(term_vars[bad]), tolower(known))]
    stop("'", paste(term_labels[bad], collapse = "', '"), "' ", if (sum(bad) > 1) "are" else "is",
         " not in the fixed effects of the model.",
         if (!anyNA(hit)) paste0(" Did you mean '", paste(hit, collapse = "', '"), "'?") else
           paste0(" Their variables are: ", paste(known, collapse = ", "), "."), call. = FALSE)
  }
  invisible(TRUE)
}

## The covariate values asked for. 'covariate' may be
##   a named list or named numeric vector, e.g. list(x = c(-1, 0, 1), z = 0.3):
##     several values give rows of K (all combinations), one value fixes the
##     covariate; a data frame gives joint points, one row of K for each of
##     its rows (per level combination);
##   an unnamed numeric vector with one value for each covariate of the
##     model, in their order;
##   covariate names, with their values in 'at' (a vector for one
##     covariate; a list, by name or in the order of 'covariate', or a matrix
##     of joint points for several); without 'at', 50 (one covariate)
##     or 30 values over the observed range.
## Returns list(fixed (named list of single values), values (named list of
## the values that define rows), varied (their names, in order), joint).
pm_cov_values <- function(covariate, at, covs, X) {
  out <- list(fixed = list(), values = list(), varied = character(0), joint = NULL)
  if (pm_is_null(covariate)) return(out)
  known <- function(nm) {
    if (any(!nm %in% covs)) {
      stop("Unknown covariate(s) in 'covariate': ", paste(setdiff(nm, covs), collapse = ", "), ". ",
           if (length(covs)) paste0("The covariates of the model are: ", paste(covs, collapse = ", "), ".") else
             "The model has no covariates.", call. = FALSE)
    }
    if (anyDuplicated(nm)) stop("A covariate is given twice in 'covariate'.", call. = FALSE)
  }
  num <- function(v, nm) {
    if (!is.numeric(v) || !length(v) || anyNA(v)) stop("The values of covariate '", nm, "' must be numbers.", call. = FALSE)
    v
  }
  if (is.data.frame(covariate) || (is.matrix(covariate) && is.numeric(covariate) && !is.null(colnames(covariate)))) {
    joint <- as.data.frame(covariate)
    known(names(joint))
    out$joint <- joint
    out$varied <- names(joint)
  } else if (is.list(covariate)) {
    if (is.null(names(covariate)) || any(!nzchar(names(covariate)))) {
      stop("A list in 'covariate' must be named by covariate, e.g. list(x = c(0, 1), z = 2).", call. = FALSE)
    }
    known(names(covariate))
    for (nm in names(covariate)) {
      v <- unique(num(covariate[[nm]], nm))
      if (length(v) == 1) out$fixed[[nm]] <- v else out$values[[nm]] <- v
    }
    out$varied <- names(out$values)
  } else if (is.numeric(covariate)) {
    if (!is.null(names(covariate)) && all(nzchar(names(covariate)))) {
      known(names(covariate))
      out$fixed <- as.list(covariate)
    } else {
      if (length(covariate) != length(covs)) {
        stop("'covariate' has ", length(covariate), " value(s) and the model has ", length(covs), " covariate(s)",
             if (length(covs)) paste0(" (", paste(covs, collapse = ", "), ")"), ". Give one value for each, in ",
             "this order, or name them: covariate = c(", if (length(covs)) covs[1] else "x", " = 1).", call. = FALSE)
      }
      out$fixed <- setNames(as.list(unname(covariate)), covs)
    }
  } else if (is.character(covariate)) {
    known(covariate)
    if (!pm_is_null(at) && (is.matrix(at) || is.data.frame(at)) && length(covariate) > 1) {
      # joint points, one row per prediction point
      joint <- as.data.frame(at)
      if (is.matrix(at) && is.null(colnames(at))) names(joint) <- covariate
      if (!all(covariate %in% names(joint))) {
        stop("The columns of 'at' must be named ", paste(covariate, collapse = ", "))
      }
      out$joint <- joint[covariate]
    } else if (!pm_is_null(at)) {
      if (!is.list(at)) at <- list(at)
      if (length(at) != length(covariate)) {
        stop("'at' must have the values of ", length(covariate), " covariate(s): ",
             paste(covariate, collapse = ", "), ".", call. = FALSE)
      }
      by_name <- !is.null(names(at)) && all(nzchar(names(at)))
      if (by_name && !setequal(names(at), covariate)) {
        stop("The names of 'at' must be those in 'covariate': ", paste(covariate, collapse = ", "), ".",
             call. = FALSE)
      }
      for (i in seq_along(covariate)) {
        v <- num(at[[if (by_name) covariate[i] else i]], covariate[i])
        out$values[[covariate[i]]] <- if (length(covariate) == 1) v else sort(unique(v))
      }
    } else if (length(covariate) == 1) {
      x <- X[[covariate]]
      out$values[[covariate]] <- sort(unique(c(seq(min(x), max(x), length = 50), x)))
    } else {
      for (i in covariate) out$values[[i]] <- sort(unique(seq(min(X[[i]]), max(X[[i]]), length = 30)))
    }
    out$varied <- covariate
  } else {
    stop("'covariate' must be numbers (values), names of covariates, or a named list of values.", call. = FALSE)
  }
  out
}

## The values of each variable in the reference grid: the levels of factors,
## the values given for covariates ('cv': pm_cov_values()), and 'cov_reduce'
## of the other covariates (their mean by default).
pm_base_values <- function(pv, cv, cov_reduce) {
  X <- pv$X
  reduce <- function(x) {
    if (identical(cov_reduce, "average")) return(mean(x, na.rm = TRUE))
    val <- cov_reduce(x)
    if (!is.numeric(val) || length(val) != 1) {
      stop("'cov_reduce' must return a single number.")
    }
    val
  }
  baselevs <- list()
  for (xname in names(X)) {
    obj <- X[[xname]]
    if (is.matrix(obj)) next                                   # matrix covariates: column means (matdat)
    baselevs[[xname]] <- if (is.factor(obj)) {
      factor(levels(obj), levels = levels(obj), ordered = is.ordered(obj))
    } else if (xname %in% pv$coerced) {
      sort(unique(obj))
    } else if (!is.null(cv$fixed[[xname]])) {
      cv$fixed[[xname]]
    } else if (!is.null(cv$values[[xname]])) {
      cv$values[[xname]]
    } else {
      reduce(obj)
    }
  }
  # Variables used only in excluded smooths: fixed at one value, never averaged.
  for (v in pv$excluded_only) baselevs[[v]] <- baselevs[[v]][1]
  baselevs
}

## Covariates that take one value within each level of a factor (with at
## least two distinct values): list(cov = list(factor =, values = named by
## level)). The first such factor is used.
pm_level_covariates <- function(X, covs, factor_names) {
  out <- list()
  for (cv in covs) {
    x <- X[[cv]]
    if (!is.numeric(x) || is.matrix(x)) next
    for (f in factor_names) {
      g <- X[[f]]
      if (!is.factor(g)) next
      ok <- !is.na(x) & !is.na(g)
      if (!any(ok)) next
      rng <- tapply(x[ok], droplevels(g[ok]), function(v) diff(range(v)))
      vals <- tapply(x[ok], droplevels(g[ok]), `[`, 1)
      if (all(rng <= 1e-10 * max(1, abs(vals))) && length(unique(vals)) > 1) {
        out[[cv]] <- list(factor = f, values = vals)
        break
      }
    }
  }
  out
}

## Set the level covariates of the rows of 'nd' to their level's value
## (levels without data keep the value they have).
pm_set_level_covariates <- function(nd, lc) {
  for (cv in names(lc)) {
    f <- lc[[cv]]$factor
    if (is.null(nd[[f]])) next
    v <- unname(lc[[cv]]$values[as.character(nd[[f]])])
    nd[[cv]] <- ifelse(is.na(v), nd[[cv]], v)
  }
  nd
}

## Means that cannot be estimated because a factor is nested in another in
## the data while the formula has them as separate effects: say so, and how
## to get the means of the combinations that exist.
pm_note_nesting <- function(undeclared, facs, weights) {
  if (!length(undeclared) || weights %in% c("cells", "flat")) return(invisible(NULL))
  said <- character(0)
  for (p in undeclared) {
    if (!any(c(p$a, p$b) %in% facs) || paste(p$b, p$a) %in% said) next     # (aliased factors: both directions)
    said <- c(said, paste(p$a, p$b))
    message("In the data '", p$b, "' is nested in '", p$a, "' (each level of ", p$b, " occurs with one level of ",
            p$a, "), and the model has them as separate effects, which cannot be told apart. To get the means, ",
            "write ", p$a, "/", p$b, " in the model formula (or use weights = \"cells\").")
  }
  invisible(NULL)
}

## Weights of the cells for the 'weights' schemes other than "equal":
## list(fun = function(grid) -> weight of each row's cell (NULL: equal), sets:
## the factors that are weighted jointly, so their combinations are needed).
pm_cell_weights <- function(weights, preddf, factor_names, nuisance) {
  if (weights == "equal" || !length(factor_names)) return(list(fun = NULL, sets = list()))
  cell_key <- function(df, vs) {
    do.call(paste, c(lapply(vs, function(v) as.character(df[[v]])), sep = "\r"))
  }
  # number of observations in each grid row's combination of the factors 'vs'
  count_of <- function(vs) {
    tab <- table(cell_key(preddf, vs))
    function(grid) {
      n <- as.vector(tab[cell_key(grid, vs)])
      n[is.na(n)] <- 0
      n
    }
  }
  switch(weights,
    # same weights for every mean: frequencies of the averaged-over cells (SAS OM)
    proportional = list(fun = if (length(nuisance)) count_of(nuisance), sets = list(nuisance)),
    # product of each averaged-over factor's own frequencies
    outer = list(fun = if (length(nuisance)) {
      cnt <- lapply(nuisance, count_of)
      function(grid) Reduce(`*`, lapply(cnt, function(f) f(grid)))
    }, sets = list()),
    # each mean weights the cells by their own frequencies
    cells = list(fun = count_of(factor_names), sets = list(factor_names)),
    # equal weights over the non-empty cells
    flat = list(fun = local({cnt <- count_of(factor_names); function(grid) as.numeric(cnt(grid) > 0)}),
                sets = list(factor_names)))
}

## The basis of the fixed effects for rows of predictor values: list(fun,
## Terms). The levels of every factor of the model frame (incl. factor(x)
## terms and character variables) are fixed, so any rows give the same columns.
pm_grid_basis <- function(info, model_data, xlev) {
  model <- info$model
  Terms <- info$Terms
  if (info$builder == "gam") {
    gam_excl <- if (length(info$gam_exclude)) info$gam_exclude else NULL
    return(list(fun = function(nd) pm_basis(model, nd, exclude = gam_excl), Terms = Terms))
  }
  contrasts <- info$contrasts
  mf_fit <- model.frame(Terms, model_data, na.action = na.omit, drop.unused.levels = TRUE)
  if (is.null(attr(Terms, "predvars"))) Terms <- attr(mf_fit, "terms")
  is_lev <- vapply(mf_fit, function(z) is.factor(z) || is.character(z), logical(1))
  mxlev <- lapply(mf_fit[is_lev], function(z) if (is.factor(z)) levels(z) else sort(unique(z)))
  mxlev <- mxlev[!names(mxlev) %in% names(xlev)]
  xlev_all <- c(xlev, mxlev)
  list(fun = function(nd) pm_basis(model, nd, terms = Terms, contrasts = contrasts, xlev = xlev_all), Terms = Terms)
}

## The sets of variables that share columns of the basis: the variables of
## each term of the model (mgcv: of each parametric term and each smooth).
pm_basis_sets <- function(info, Terms) {
  if (info$builder != "gam") return(pm_term_sets(Terms))
  sm <- info$model$smooth
  c(pm_term_sets(stats::delete.response(info$model$pterms)),
    lapply(sm, function(s) setdiff(c(s$term, s$by), "NA")))
}

## The position, in each dimension of the grid, of the values of the first
## row of the data: an observed combination, which every restriction on the
## cells (observed cells, nested factors, frequencies) accepts.
pm_ref_positions <- function(dims, X) {
  lapply(dims, function(d) {
    if (length(d$vars) != 1 || nrow(d$values) == 1 || !d$vars %in% names(X)) return(1L)
    k <- match(as.character(X[[d$vars]][1]), as.character(d$values[[1]]))
    if (is.na(k)) 1L else k
  })
}

## Null space of a rank-deficient design X, for the estimability checks of
## Kmatrix(), the empty-cell contrasts, the pairs of predictmeans() and
## pm_slopes(). Computed on X with its columns scaled to unit length, so that
## the check does not depend on the units of a covariate (a site mean in
## metres or in millimetres). Returns D^-1 N (exact null vectors of X; D the
## column lengths) with attribute "scale" = D, or NULL for a full-rank X.
pm_null_space <- function(X) {
  X <- as.matrix(X)
  d <- sqrt(colSums(X^2))
  d[!is.finite(d) | d == 0] <- 1
  ev <- eigen(crossprod(sweep(X, 2, d, "/")), symmetric = TRUE)
  tol <- max(dim(X)) * max(ev$values) * .Machine$double.eps * 100
  N <- ev$vectors[, ev$values <= tol, drop = FALSE]
  if (!ncol(N)) return(NULL)
  N <- N / d
  rownames(N) <- colnames(X)
  structure(N, scale = stats::setNames(d, colnames(X)))
}

## Rows of L (linear functions of the coefficients) that are estimable: their
## projection 'np' = L N on the null space (pm_null_space(), with its
## "scale") is zero relative to L in the scaled columns (L D^-1).
pm_estimable_rows <- function(np, L) {
  d <- attr(np, "scale")
  L <- as.matrix(L)
  # columns of L that are not in the design (none, or coefficients already
  # dropped) count as 0
  cn <- if (is.null(d)) NULL else intersect(names(d), colnames(L))
  Ls <- if (is.null(d)) L else sweep(L[, cn, drop = FALSE], 2, d[cn], "/")
  apply(abs(as.matrix(np)), 1, max) <= 1e-6 * apply(abs(Ls), 1, max)
}

