Kmatrix <- function(model, modelterm, covariate=NULL, covariateV=NULL, data=NULL,
                    prtnum=FALSE, cov_reduce=mean, exclude=NULL,
                    weights=c("equal", "proportional", "outer", "cells", "flat")) {

  weights <- match.arg(weights)
  info <- pm_model_info(model, data, exclude)
  model <- info$model
  Terms <- info$Terms
  yname <- info$yname

  rec <- pm_recover_data(info)
  model_data <- rec$data
  nm <- rec$needed

  # Raw predictor values (one column per variable) for the rows used in the fit.
  form <- as.formula(paste("~", paste(pm_bq(nm), collapse = "+")), env = info$envir)
  X <- model.frame(form, model_data, na.action = na.omit, drop.unused.levels = TRUE)
  for (v in names(X)) {
    if (is.character(X[[v]]) || is.logical(X[[v]])) X[[v]] <- factor(X[[v]])
  }
  preddf <- X

  # Variables wrapped in factor() etc. in the formula are treated as factors
  # with their observed values as levels.
  coerced <- intersect(pm_coerced_vars(formula(Terms)), names(X))
  excluded_only <- intersect(info$excluded_only, names(X))

  term_labels <- unlist(strsplit(modelterm, "\\:"))
  term_vars <- pm_term_vars(term_labels)

  cov_given <- !pm_is_null(covariate)
  if (!cov_given || !setequal(covariate, modelterm)) {
    is_fac <- vapply(term_vars, function(v) {
      v %in% names(X) && (is.factor(X[[v]]) || v %in% coerced)
    }, logical(1))
    if (any(!is_fac)) {
      stop(paste(term_labels[!is_fac], "must be factor(s) or not in the fixed effects!"))
    }
    if (any(term_vars %in% excluded_only)) {
      stop(paste("Factor(s)", paste(term_labels[term_vars %in% excluded_only], collapse = ", "),
                 "only appear in smooth terms excluded from prediction (see 'exclude')."))
    }
  }

  average_cov <- identical(cov_reduce, "average")
  if (!average_cov && !is.function(cov_reduce)) {
    stop("'cov_reduce' must be a function such as mean or median, or \"average\".")
  }
  reduce <- function(x) {
    if (average_cov) return(mean(x, na.rm = TRUE))
    val <- cov_reduce(x)
    if (!is.numeric(val) || length(val) != 1) {
      stop("'cov_reduce' must return a single number.")
    }
    val
  }

  baselevs <- xlev <- matdat <- list()
  for (xname in names(X)) {
    obj <- X[[xname]]
    if (is.factor(obj)) {
      xlev[[xname]] <- levels(obj)
      baselevs[[xname]] <- factor(levels(obj), levels = levels(obj), ordered = is.ordered(obj))
    } else if (is.matrix(obj)) {
      # Matrix covariates: use column means
      matdat[[xname]] <- apply(obj, 2, mean, na.rm = TRUE)
    } else if (xname %in% coerced) {
      baselevs[[xname]] <- sort(unique(obj))
    } else {
      baselevs[[xname]] <- reduce(obj)
    }
  }
  # Variables used only in excluded smooths: fixed at one value, never averaged.
  for (v in excluded_only) baselevs[[v]] <- baselevs[[v]][1]

  factor_names <- setdiff(c(names(xlev), coerced), excluded_only)
  covlevname <- setdiff(names(baselevs), c(names(xlev), coerced))
  covlevname <- setdiff(covlevname, c(info$offset_vars, excluded_only))

  joint <- NULL
  if (cov_given && is.numeric(covariate)) {
    baselevs[covlevname] <- as.list(covariate)
  }
  if (cov_given && is.character(covariate) && all(covariate %in% covlevname)) {
    if (length(covariate) == 1) {
      if (!pm_is_null(covariateV) && is.vector(covariateV)) {
        baselevs[[covariate]] <- covariateV
      } else {
        baselevs[[covariate]] <- sort(unique(c(seq(min(X[[covariate]]), max(X[[covariate]]), length=50),
                                               X[[covariate]])))
      }
    } else {
      if (!pm_is_null(covariateV) && is.matrix(covariateV)) {
        # A matrix gives joint points, one row per prediction point.
        joint <- as.data.frame(covariateV)
        if (is.null(colnames(covariateV))) names(joint) <- covariate
        if (!all(covariate %in% names(joint))) {
          stop("The columns of 'covariateV' must be named ", paste(covariate, collapse = ", "))
        }
        joint <- joint[covariate]
      } else if (!pm_is_null(covariateV) && is.list(covariateV)) {
        for (i in seq_along(covariate))
          baselevs[[covariate[i]]] <- sort(unique(covariateV[[i]]))
      } else {
        for (i in covariate)
          baselevs[[i]] <- sort(unique(c(seq(min(X[[i]]), max(X[[i]]), length=30))))
      }
    }
  }

  # These are the variables defining the rows of K.
  facs <- term_vars
  if (cov_given && all(is.character(covariate), all(!covariate %in% facs))) {
    facs <- c(facs, covariate)
  }
  if (any(!facs %in% names(baselevs))) {
    stop(paste("Unknown factor(s) in specification:", paste(facs[!facs %in% names(baselevs)], collapse=" ")))
  }

  # Optionally average over the observed values of the covariates that are
  # not varied, instead of fixing them at one value.
  profiles <- prof_w <- NULL
  avg_vars <- character(0)
  if (average_cov) {
    avg_vars <- setdiff(covlevname, facs)
    if (length(avg_vars)) {
      pdat <- X[avg_vars]
      pkey <- do.call(paste, c(lapply(pdat, function(z) format(z, digits = 15)), sep = "\r"))
      first <- !duplicated(pkey)
      profiles <- pdat[first, , drop = FALSE]
      rownames(profiles) <- NULL
      prof_w <- as.vector(table(factor(pkey, levels = pkey[first]))) / length(pkey)
    }
  }

  if (prtnum && length(covlevname) != 0) {
    fixed_cov <- setdiff(covlevname, c(facs, avg_vars))
    if (length(fixed_cov)) {
      message("\n The predicted means are estimated at \n")
      message(paste(utils::capture.output(print(round(unlist(baselevs[fixed_cov]), 4))), collapse = "\n"))
    }
    if (length(avg_vars)) {
      message("\n The predicted means are averaged over the observed values of ",
              paste(avg_vars, collapse = ", "))
    }
  }

  # Weights for averaging over the factors that are not in modelterm.
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
  nuisance <- setdiff(factor_names, facs)
  present <- NULL
  if (weights != "equal" && length(factor_names) != 0) {
    present <- switch(weights,
      # same weights for every mean: frequencies of the averaged-over cells (SAS OM)
      proportional = if (length(nuisance)) count_of(nuisance),
      # product of each averaged-over factor's own frequencies
      outer = if (length(nuisance)) function(grid) Reduce(`*`, lapply(nuisance, function(v) count_of(v)(grid))),
      # each mean weights the cells by their own frequencies
      cells = count_of(factor_names),
      # equal weights over the non-empty cells
      flat = function(grid) as.numeric(count_of(factor_names)(grid) > 0))
  }

  # Missing factor combinations: only average over the observed cells.
  # Only the factors that interact matter: blocks or other main-effect-only
  # factors would make every cell look single (one observation per cell),
  # hiding structurally empty combinations such as nested treatments.
  if (weights == "equal" && length(factor_names) != 0) {
    fac_mat <- attr(Terms, "factors")
    ord <- attr(Terms, "order")
    inter_rows <- if (length(fac_mat) && any(ord > 1))
      rownames(fac_mat)[rowSums(fac_mat[, ord > 1, drop = FALSE] > 0) > 0] else character(0)
    inter_vars <- factor_names[factor_names %in% inter_rows |
                                 vapply(factor_names, function(v) any(grepl(paste0("\\b", v, "\\b"), inter_rows)), logical(1))]
    if (length(inter_vars) < 2) inter_vars <- factor_names
    ndf <- data.frame(table(preddf[, inter_vars, drop = FALSE]))
    max_ord <- max(c(0, ord))
    if ((!cov_given || !setequal(covariate, modelterm)) && max_ord > 1 &&
        length(inter_vars) > 1 && any(ndf$Freq == 0) &&
        !setequal(unique(ndf$Freq), c(0, 1))) {
      warning("Missing treatments' combination appeared, predicted means maybe misleading!")
      seen <- unique(cell_key(preddf, inter_vars))
      present <- function(grid) cell_key(grid, inter_vars) %in% seen
    }
  }

  rg <- pm_make_refgrid(facs, baselevs, matdat, joint, profiles, prof_w, present)

  # Basis evaluated row by row on the grid, through pm_basis() (see ?pm_coef).
  if (info$builder == "gam") {
    gam_excl <- if (length(info$gam_exclude)) info$gam_exclude else NULL
    basis <- function(nd) pm_basis(model, nd, exclude = gam_excl)
  } else {
    contrasts <- info$contrasts
    # Levels of every factor in the model frame (incl. factor(x) terms and
    # character variables), so any subset of grid rows gives the same columns.
    mf_fit <- model.frame(Terms, model_data, na.action = na.omit, drop.unused.levels = TRUE)
    if (is.null(attr(Terms, "predvars"))) Terms <- attr(mf_fit, "terms")
    is_lev <- vapply(mf_fit, function(z) is.factor(z) || is.character(z), logical(1))
    mxlev <- lapply(mf_fit[is_lev], function(z) if (is.factor(z)) levels(z) else sort(unique(z)))
    mxlev <- mxlev[!names(mxlev) %in% names(xlev)]
    xlev_all <- c(xlev, mxlev)
    basis <- function(nd) pm_basis(model, nd, terms = Terms, contrasts = contrasts, xlev = xlev_all)
  }

  K <- as.matrix(pm_refgrid_average(rg, basis))

  fctnames <- rg$target
  names(fctnames)[seq_along(term_vars)] <- term_labels
  rnK <- do.call("paste", c(fctnames[, seq_along(facs), drop = FALSE], sep=":"))
  rownames(K) <- rnK

  if (length(setdiff(info$coef_names, colnames(K))) != 0) {
    stop(paste("You may need drop empty levels in factor", modelterm))
  }

  # Estimability: a row of K is estimable if it is orthogonal to the null
  # space of the design matrix used in the fit. Only rank-deficient fits
  # (aliased or dropped coefficients) have a null space, so the check costs
  # nothing for full-rank models.
  estimable <- rep(TRUE, nrow(K))
  if (info$builder != "gam" && info$rank_deficient(colnames(K))) {
    Xfit <- basis(model_data)
    XtX <- crossprod(Xfit)
    ev <- eigen(XtX, symmetric = TRUE)
    tol <- max(dim(Xfit)) * max(ev$values) * .Machine$double.eps * 100
    N <- ev$vectors[, ev$values <= tol, drop = FALSE]
    if (ncol(N)) {
      Kfull <- K[, colnames(Xfit), drop = FALSE]
      nonest <- abs(Kfull %*% N)
      estimable <- apply(nonest, 1, max) <= 1e-6 * pmax(1, apply(abs(Kfull), 1, max))
    }
  }
  estimable[is.na(estimable) | rg$empty] <- FALSE
  names(estimable) <- rnK

  K <- K[, info$coef_names, drop = FALSE]
  attr(K, "estimable") <- estimable

  # Covariates that are fixed (not varied) are reported with the value used.
  cov_report <- setdiff(covlevname, names(fctnames))
  if (length(cov_report)) {
    cov_vals <- as.matrix(rg$W %*% as.matrix(rg$grid[cov_report]))
    if (any(rg$empty)) cov_vals[rg$empty, ] <- NaN
    colnames(cov_vals) <- cov_report
    fctnames <- cbind(fctnames, cov_vals)
  }
  rownames(fctnames) <- NULL

  return(list(K=K, fctnames=fctnames, response=yname, preddf=preddf, refgrid=rg, data=model_data))
}
