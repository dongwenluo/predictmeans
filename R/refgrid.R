## Reference-grid engine used by Kmatrix() and semipred().
##
## A reference grid is a data frame of predictor values ("grid") plus a sparse
## weight matrix W. Row t of W averages the grid rows that belong to target
## level combination t (equal weights over the other factors' cells, and
## optional frequency weights over observed covariate values). Any basis that
## can be evaluated row by row on new data -- model.matrix(), an mgcv
## lpmatrix, or the smooth Z matrices of semireg() -- is averaged the same way
## with pm_refgrid_average(), so fixed and smooth parts always line up.

## TRUE for NULL, zero-length input, or the legacy "NULL"/"" strings.
pm_is_null <- function(x) {
  is.null(x) || length(x) == 0 ||
    (is.character(x) && all(x %in% c("NULL", "")))
}

## Variables wrapped in factor(), as.factor(), ordered() or as.ordered()
## anywhere in an expression, e.g. y ~ factor(dose) + x gives "dose".
pm_coerced_vars <- function(expr) {
  out <- character(0)
  walk <- function(e) {
    if (is.call(e)) {
      fn <- e[[1]]
      if (is.name(fn) &&
          as.character(fn) %in% c("factor", "as.factor", "ordered", "as.ordered")) {
        out <<- c(out, all.vars(e))
      } else {
        for (a in as.list(e)[-1]) walk(a)
      }
    }
  }
  walk(expr)
  unique(out)
}

## Map each component of a term label such as "factor(dose):A" to the
## variable it refers to ("dose", "A").
pm_term_vars <- function(labels) {
  vapply(labels, function(lb) {
    e <- tryCatch(str2lang(lb), error = function(err) NULL)
    if (is.call(e)) {
      v <- all.vars(e)
      if (length(v) == 1) return(v)
    }
    lb
  }, character(1), USE.NAMES = FALSE)
}

## Backtick-quote names so they can be pasted into a formula.
pm_bq <- function(x) paste0("`", gsub("`", "\\\\`", x), "`")

## Class-specific facts needed to build a reference grid.
pm_model_info <- function(model, data = NULL, exclude = NULL) {
  extra_vars <- character(0)
  data_first <- FALSE
  if (inherits(model, "semireg")) {
    # Smooth-only covariates live in the smoothZ list, not the fixed formula.
    if (is.null(data)) data <- model$data
    data_first <- TRUE
    extra_vars <- setdiff(model$sm_vars, model$fomul_vars)
    model <- model$semer
  }

  builder <- "mm"
  if (inherits(model, "gam")) {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is needed for 'gam' models.")
    }
    builder <- "gam"
    thecall <- model$call
    contrasts <- NULL
    coef_names <- names(coef(model))
  } else if (inherits(model, "mer") || inherits(model, "merMod")) {
    if (!lme4::isLMM(model) && !lme4::isGLMM(model)) {
      stop("Can't handle a nonlinear mixed model")
    }
    thecall <- slot(model, "call")
    mm <- model.matrix(model)
    contrasts <- attr(mm, "contrasts")
    coef_names <- colnames(mm)
  } else if (inherits(model, "lme") || inherits(model, "gls")) {
    thecall <- model$call
    contrasts <- model$contrasts
    coef_names <- colnames(model.matrix(model))
  } else if (inherits(model, "lm") || inherits(model, "glmmTMB")) {
    thecall <- model$call
    mm <- model.matrix(model)
    contrasts <- attr(mm, "contrasts")
    coef_names <- colnames(mm)
  } else {
    # Any other model class, through the adapter generics (see ?pm_coef).
    builder <- "custom"
    coef_try <- tryCatch(pm_coef(model), error = function(e) NULL)
    if (!is.numeric(coef_try) || is.null(names(coef_try))) {
      stop("Can't handle a model of class ", sQuote(class(model)[1]),
           ": it needs named fixed-effect estimates from coef() or a pm_coef() method (see ?pm_coef).",
           call. = FALSE)
    }
    thecall <- tryCatch(stats::getCall(model), error = function(e) NULL)
    contrasts <- pm_contrasts(model)
    coef_names <- names(coef_try)
  }

  # A fit is rank deficient if coefficients are aliased (NA) or design
  # columns were dropped before fitting (lme4, glmmTMB).
  coef_vals <- tryCatch(pm_coef(model), error = function(e) NULL)
  coef_na <- any(is.na(coef_vals))
  rank_deficient <- function(full_cols) {
    builder != "gam" && (coef_na || length(setdiff(full_cols, coef_names)) > 0)
  }

  Terms_full <- terms(model)
  envir <- attr(Terms_full, ".Environment")
  if (is.null(envir)) envir <- environment(formula(model))
  vars_expr <- attr(Terms_full, "variables")

  yname <- NULL
  if (attr(Terms_full, "response") > 0) {
    yname <- as.character(vars_expr)[[attr(Terms_full, "response") + 1]]
  }

  offset_vars <- character(0)
  if (!is.null(attr(Terms_full, "offset"))) {
    offset_vars <- unique(unlist(lapply(attr(Terms_full, "offset"),
                                        function(i) all.vars(vars_expr[[i + 1]]))))
  }

  Terms <- delete.response(Terms_full)
  rhs_vars <- all.vars(formula(Terms))
  resp_vars <- setdiff(all.vars(formula(Terms_full)), rhs_vars)

  # mgcv: smooth terms that are excluded from predictions (random effects by
  # default), and variables that only appear in those terms.
  gam_exclude <- character(0)
  excluded_only <- character(0)
  if (builder == "gam") {
    sm <- model$smooth
    sm_labels <- vapply(sm, function(s) s$label, character(1))
    if (is.null(exclude)) {
      is_random <- vapply(sm, function(s) {
        inherits(s, "random.effect") || inherits(s, "fs.interaction")
      }, logical(1))
      gam_exclude <- sm_labels[is_random]
    } else {
      bad <- setdiff(exclude, sm_labels)
      if (length(bad)) {
        stop("Unknown smooth term(s) in 'exclude': ", paste(bad, collapse = ", "),
             ". Available: ", paste(sm_labels, collapse = ", "))
      }
      gam_exclude <- exclude
    }
    sm_vars <- function(s) setdiff(c(s$term, s$by), "NA")
    kept_vars <- unique(c(all.vars(delete.response(model$pterms)),
                          unlist(lapply(sm[!sm_labels %in% gam_exclude], sm_vars))))
    excl_vars <- unique(unlist(lapply(sm[sm_labels %in% gam_exclude], sm_vars)))
    excluded_only <- setdiff(excl_vars, kept_vars)
  }

  list(model = model, builder = builder, thecall = thecall, contrasts = contrasts,
       coef_names = coef_names, rank_deficient = rank_deficient,
       Terms = Terms, envir = envir, yname = yname,
       offset_vars = offset_vars, rhs_vars = rhs_vars, resp_vars = resp_vars,
       extra_vars = extra_vars, data = data, data_first = data_first,
       gam_exclude = gam_exclude, excluded_only = excluded_only)
}

## Recover the raw (untransformed) predictor data for the rows used in the fit.
pm_recover_data <- function(info) {
  model <- info$model
  needed <- unique(c(info$rhs_vars, info$extra_vars))
  wanted <- unique(c(needed, info$resp_vars))
  safe <- function(expr) tryCatch(expr, error = function(e) NULL)
  frame <- safe(model.frame(model))
  is_nlme <- inherits(model, "lme") || inherits(model, "gls")
  is_mer <- inherits(model, "mer") || inherits(model, "merMod")

  # Possible sources of the raw data, tried in order (evaluated lazily).
  cands <- list()
  if (!is.null(info$data) && (info$data_first || !(is_nlme || is_mer))) {
    cands$data <- function() info$data
  }
  cands$pm_data <- function() pm_data(model)
  if (is_nlme) cands$getData <- function() pm_get_data(model)
  if (is_mer) cands$getData <- function() nlme::getData(model)
  if (!is.null(info$data) && is.null(cands$data)) cands$data <- function() info$data
  if (is.null(info$data)) {
    if (!is.null(info$thecall$data)) {
      cands$call <- function() eval(info$thecall$data, info$envir)
    }
    cands$frame <- function() frame
    cands$env <- function() {
      get_all_vars(as.formula(paste("~", paste(pm_bq(needed), collapse = "+")), env = info$envir))
    }
  }

  # Names that are not data columns but scalars in the formula environment
  # (e.g. 'pi' in I(pi * x)) are constants, not variables.
  is_const <- function(v) {
    val <- safe(get(v, envir = info$envir))
    !is.null(val) && length(val) == 1 && !is.function(val)
  }

  dat <- NULL
  seen <- character(0)
  constants <- character(0)
  for (src in names(cands)) {
    d <- safe(cands[[src]]())
    if (!(is.data.frame(d) || is.list(d))) next
    if (!is.data.frame(d)) d <- as.data.frame(d, check.names = FALSE)
    seen <- union(seen, names(d))
    miss <- setdiff(needed, names(d))
    if (length(miss) && all(vapply(miss, is_const, logical(1)))) constants <- miss
    if (all(setdiff(needed, constants) %in% names(d))) {
      dat <- d
      from_getData <- src == "getData"
      break
    }
  }
  if (is.null(dat)) {
    missing_v <- setdiff(needed, seen)
    stop("Can't find the variable(s) ", paste(sQuote(missing_v), collapse = ", "),
         " used in the model. Please supply the data set used in model fitting via 'data'.")
  }
  needed <- setdiff(needed, constants)
  wanted <- setdiff(wanted, constants)

  # Keep only the rows that were used in the fit.
  if (!is.null(frame) && !is_nlme && nrow(frame) <= nrow(dat) &&
      all(rownames(frame) %in% rownames(dat))) {
    dat <- dat[rownames(frame), , drop = FALSE]
  } else if (!is.null(info$thecall$subset) && !(from_getData && is_nlme)) {
    keep <- safe(eval(info$thecall$subset, dat, info$envir))
    if (is.logical(keep) && length(keep) == nrow(dat)) {
      dat <- dat[keep & !is.na(keep), , drop = FALSE]
    }
  }

  dat <- dat[, intersect(wanted, names(dat)), drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  list(data = droplevels(dat), needed = needed, constants = constants)
}

## Build the reference grid.
##   facs        variables defining the target rows (factors and/or varied covariates)
##   baselevs    named list of values for each single variable dimension
##   joint       optional data frame of jointly specified covariate points
##   profiles    optional data frame of observed covariate profiles to average over
##   prof_w      weights for 'profiles' (sum to 1)
##   present     function(grid) -> weight for each grid row's factor cell (logical
##               TRUE/FALSE for usable/unusable cells, or numeric cell weights)
pm_make_refgrid <- function(facs, baselevs, matdat, joint = NULL,
                            profiles = NULL, prof_w = NULL, present = NULL) {
  dims <- list()
  for (v in names(baselevs)) {
    dims[[v]] <- list(vars = v, values = setNames(data.frame(baselevs[[v]],
                                                             stringsAsFactors = FALSE), v))
  }
  if (!is.null(joint)) {
    for (v in names(joint)) dims[[v]] <- NULL
    dims[[".joint"]] <- list(vars = names(joint), values = joint)
  }
  if (!is.null(profiles)) {
    for (v in names(profiles)) dims[[v]] <- NULL
    dims[[".profile"]] <- list(vars = names(profiles), values = profiles, w = prof_w)
  }

  sizes <- vapply(dims, function(d) nrow(d$values), integer(1))
  idx <- do.call(expand.grid, c(lapply(sizes, seq_len), KEEP.OUT.ATTRS = FALSE))
  grid <- do.call(cbind, unname(lapply(names(dims), function(nm) {
    dims[[nm]]$values[idx[[nm]], , drop = FALSE]
  })))
  rownames(grid) <- NULL
  for (nm in names(matdat)) {
    grid[[nm]] <- matrix(rep(matdat[[nm]], each = nrow(grid)), nrow = nrow(grid),
                         dimnames = list(NULL, names(matdat[[nm]])))
  }

  # Target dimensions: first one varies slowest, as in Kmatrix() 1.x.
  tdims <- names(dims)[vapply(dims, function(d) any(d$vars %in% facs), logical(1))]
  tdims <- tdims[order(vapply(tdims, function(nm) min(match(dims[[nm]]$vars, facs), na.rm = TRUE),
                              numeric(1)))]
  tsizes <- sizes[tdims]
  strides <- rev(cumprod(rev(c(tsizes[-1], 1))))
  key <- rep(1L, nrow(grid))
  for (k in seq_along(tdims)) key <- key + (idx[[tdims[k]]] - 1L) * as.integer(strides[k])
  tidx <- do.call(expand.grid, c(lapply(rev(tsizes), seq_len), KEEP.OUT.ATTRS = FALSE))
  tidx <- tidx[, rev(seq_along(tsizes)), drop = FALSE]
  names(tidx) <- tdims
  target <- do.call(cbind, unname(lapply(tdims, function(nm) {
    dims[[nm]]$values[tidx[[nm]], , drop = FALSE]
  })))
  target <- target[, intersect(facs, names(target)), drop = FALSE]
  rownames(target) <- NULL

  w <- rep(1, nrow(grid))
  if (!is.null(profiles)) w <- w * prof_w[idx[[".profile"]]]
  if (!is.null(present)) w <- w * as.numeric(present(grid))
  tot <- as.vector(tapply(w, factor(key, levels = seq_len(nrow(target))), sum))
  tot[is.na(tot)] <- 0
  x <- as.vector(w / tot[key])
  keep <- w > 0
  W <- Matrix::sparseMatrix(i = key[keep], j = which(keep), x = x[keep],
                            dims = c(nrow(target), nrow(grid)))
  structure(list(grid = grid, W = W, target = target, empty = as.vector(tot == 0)),
            class = "pm_refgrid")
}

## Average a row-wise basis over the reference grid: W %*% FUN(grid).
## Large grids are processed in chunks of rows to bound memory.
pm_refgrid_average <- function(rg, FUN, chunk = 20000L) {
  n <- nrow(rg$grid)
  starts <- seq(1L, n, by = chunk)
  out <- NULL
  for (s in starts) {
    rows <- s:min(n, s + chunk - 1L)
    B <- FUN(rg$grid[rows, , drop = FALSE])
    part <- rg$W[, rows, drop = FALSE] %*% B
    out <- if (is.null(out)) part else out + part
  }
  if (any(rg$empty)) out[rg$empty, ] <- NaN
  colnames(out) <- colnames(B)
  out
}

## Stop with a clear message when an optional (Suggests) package is missing.
pm_need <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is needed for %s. Please install it.", pkg, what), call. = FALSE)
  }
  invisible(TRUE)
}

## Which points lie inside a polygon (two-column matrices or data frames).
## Same winding-angle method as HRW::pointsInPoly() (GPL >= 2, Wand et al.),
## but an open polygon is closed by joining its last vertex to its first.
pm_points_in_poly <- function(points, polygon) {
  points <- as.matrix(points)
  polygon <- as.matrix(polygon)
  nv <- nrow(polygon)
  if (nv > 1 && all(polygon[1, ] == polygon[nv, ])) {
    polygon <- polygon[-nv, , drop = FALSE]
    nv <- nv - 1
  }
  x <- -outer(points[, 1], polygon[, 1], "-")
  y <- -outer(points[, 2], polygon[, 2], "-")
  i2 <- c(seq_len(nv)[-1], 1)
  theta <- atan2(x * y[, i2, drop = FALSE] - y * x[, i2, drop = FALSE],
                 x * x[, i2, drop = FALSE] + y * y[, i2, drop = FALSE])
  theta <- abs(rowSums(theta))
  as.vector(abs(theta - 2 * pi) < 1e-06)
}
