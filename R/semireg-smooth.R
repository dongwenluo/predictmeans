## Smooth terms written in the formula of semireg() / semireg_tmb()
## (dev/review-semireg.md, step 3):
##   y ~ A + s(x, by = A) + s(x1, x2) + (1 | Block)
##
## The smooths are set up by mgcv (gam(fit = FALSE): bases, knots and the
## identifiability constraints, as for a gam) and turned into fixed and
## random parts with mgcv::smooth2random() (as gamm4 does). Each random part
## is a block of iid random effects: one placeholder (1 | name) in the lme4
## or glmmTMB formula, its Z replaced by the block (R/semireg-engine.R). The
## fixed part (the null space of the penalty: x, f:x, ...) is added to the
## formula as ordinary terms, so that Kmatrix() and predictmeans() treat it
## as any fixed effect; terms the user already wrote are not doubled.
## Prediction: PredictMat() on new data times the same linear map that gave
## each block on the data.

## The smooth terms of a formula: list(formula (fixed part, with the lme4
## random terms and without the smooths), groups (character vectors of
## smooth calls, each group set up by one gam(fit = FALSE) call), fixed
## (fixed terms of expanded smooths), derived (by-factors made from several
## factors), vars (variables of the smooths, never centred)) or NULL when
## there is none. 'data' tells factors from covariates.
##
## A factor (or a factor term) among the arguments of s() is expanded as an
## ANOVA-like split (dev/review-semireg.md, decision 2):
##   s(x, f)       common curve s(x) + a deviation curve per level of f, one
##                 variance for all of them
##   s(x, f1:f2)   the same with the cells of f1:f2
##   s(x, f1*f2)   common curve + deviations for f1, f2 and f1:f2, a
##                 variance each
##   s(x1, x2, f)  the same with surfaces
## The fixed part: x + f + f:x (with the factor terms of the formula: x +
## f1*f2 + f1:x + f2:x + f1:f2:x). Each piece is set up on its own, so no
## constraint is put between the common curve and the deviations: the mixed
## model separates them by their variances, as nested random effects.
pm_sm_parse <- function(formula, data = NULL) {
  # glmmTMB covariance-structure terms (ar1(), ou(), ...) are kept as written
  specials <- pm_special_terms(formula)
  bars <- pm_plain_bars(formula, specials)
  deparse1l <- function(e) paste(deparse(e, width.cutoff = 500L), collapse = " ")
  # random curves (s(x) | g): (x | g), the straight lines per level with
  # correlated intercepts and slopes, plus the wiggly part of the smooth per
  # level with one variance (Durban et al., 2005; HRW section 9.2)
  rcurves <- list()
  for (i in seq_along(bars)) {
    lhs <- bars[[i]][[2]]
    if (!(is.call(lhs) && identical(as.character(lhs[[1]]), "s"))) next
    if (!identical(as.character(bars[[i]][[1]]), "|")) stop("Write a random curve with '|': (s(x) | g).", call. = FALSE)
    a <- as.list(lhs)[-1]
    an <- if (is.null(names(a))) rep("", length(a)) else names(a)
    if ("by" %in% an) stop("A random curve (s(x) | g) cannot have 'by'.", call. = FALSE)
    covs <- vapply(a[an == ""], deparse1l, "")
    if (is.null(lhs$bs)) lhs$bs <- "os"
    rcurves[[length(rcurves) + 1L]] <- list(call = deparse1l(lhs), group = deparse1l(bars[[i]][[3]]),
                                           vars = c(all.vars(lhs), all.vars(bars[[i]][[3]])))
    bars[[i]] <- call("|", str2lang(paste(covs, collapse = " + ")), bars[[i]][[3]])
  }
  f_fixed <- pm_nobars(formula)
  is_sm <- function(e) is.call(e) && as.character(e[[1]])[1] %in% c("s", "te", "t2", "ti")
  rhs_terms <- attr(stats::terms(f_fixed, specials = c("s", "te", "t2", "ti")), "variables")[-1]
  sm_calls <- Filter(is_sm, as.list(rhs_terms))
  if (!length(sm_calls) && !length(rcurves)) return(NULL)
  # a smooth may only enter as a main term (not s(x):A)
  labs <- attr(stats::terms(f_fixed), "term.labels")
  sm_lab <- vapply(sm_calls, function(e) paste(deparse(e, width.cutoff = 500L), collapse = " "), "")
  bad <- labs[!labs %in% sm_lab & vapply(labs, function(l) any(vapply(sm_lab, grepl, logical(1), x = l,
                                                                         fixed = TRUE)), logical(1))]
  if (length(bad)) stop("A smooth can only be a term of its own: ", paste(bad, collapse = ", "), call. = FALSE)
  is_factor <- function(e) {
    if (is.call(e)) return(as.character(e[[1]]) %in% c(":", "*"))
    is.name(e) && !is.null(data) && as.character(e) %in% names(data) &&
      (is.factor(data[[as.character(e)]]) || is.character(data[[as.character(e)]]) ||
         is.logical(data[[as.character(e)]]))
  }
  plain <- character(0)
  groups <- list()
  fixed <- character(0)
  derived <- list()
  vars <- character(0)
  taken <- names(data)
  used <- character(0)
  for (e in sm_calls) {
    # tensor products: mgcv's pre-fit design of t2() is not in the span of its
    # prediction matrix (mgcv 1.9-1), so its split into fixed and random parts
    # cannot be checked yet (dev/review-semireg.md); te() and ti() cannot be
    # written as independent random effects at all
    if (as.character(e[[1]]) %in% c("te", "t2", "ti")) {
      stop("Tensor-product smooths (te(), t2(), ti()) are not available in semireg() yet; use s(x1, x2) ",
           "(isotropic thin plate) or s(x1, x2, bs = \"ztps\").", call. = FALSE)
    }
    args <- as.list(e)[-1]
    nms <- if (is.null(names(args))) rep("", length(args)) else names(args)
    bs <- if ("bs" %in% nms) eval(args[["bs"]]) else "tp"
    unnamed <- args[nms == ""]
    fac <- if (any(bs %in% c("fs", "re", "sz", "mrf"))) logical(length(unnamed)) else
      vapply(unnamed, is_factor, logical(1))
    vars <- c(vars, unlist(lapply(unnamed, all.vars)), if ("by" %in% nms) all.vars(args[["by"]]))
    if (!any(fac)) {
      plain <- c(plain, deparse1l(e))
      next
    }
    if (sum(fac) > 1) stop("Give the factors of a smooth as one term, e.g. s(x, f1*f2) or s(x, f1:f2): ",
                           deparse1l(e), call. = FALSE)
    if ("by" %in% nms) stop("A smooth with a factor term cannot also have 'by': ", deparse1l(e), call. = FALSE)
    covs <- vapply(unnamed[!fac], deparse1l, "")
    if (!length(covs)) stop("A smooth needs a covariate besides its factor term: ", deparse1l(e), call. = FALSE)
    fterm <- unnamed[fac][[1]]
    tl <- attr(stats::terms(stats::as.formula(paste("~", deparse1l(fterm)))), "term.labels")
    other <- args[nms != "" & nms != "id"]
    opt <- if (length(other)) paste0(", ", paste(names(other), "=", vapply(other, deparse1l, ""), collapse = ", ")) else ""
    cov_txt <- paste(covs, collapse = ", ")
    # each piece gets its own id "pm:<block name>" (unique): the pieces are
    # found by it, and it names their block
    uniq <- function(nm) {
      while (nm %in% used) nm <- paste0(nm, "_2")
      used <<- c(used, nm)
      nm
    }
    common <- paste0("s(", cov_txt, opt, ", id = \"pm:", uniq(pm_sm_name(paste0("s(", cov_txt, ")"), character(0))), "\")")
    groups <- c(groups, list(common))
    fixed <- c(fixed, covs)
    for (t in tl) {
      fs <- strsplit(t, ":", fixed = TRUE)[[1]]
      byv <- if (length(fs) == 1) fs else {
        nm <- paste(fs, collapse = "_")
        while (nm %in% taken) nm <- paste0(nm, "_by")
        taken <- c(taken, nm)
        derived[[nm]] <- fs
        nm
      }
      blk <- uniq(pm_sm_name(paste0("s(", cov_txt, "|", t, ")"), character(0)))
      groups <- c(groups, list(paste0("s(", cov_txt, ", by = ", byv, opt, ", id = \"pm:", blk, "\")")))
      fixed <- c(fixed, t, paste0(t, ":", covs))
    }
  }
  if (length(plain)) groups <- c(list(plain), groups)
  ig <- mgcv::interpret.gam(f_fixed)
  pf <- ig$pf
  rhs <- c(paste(deparse(pf[[3]], width.cutoff = 500L), collapse = " "),
           vapply(bars, function(b) paste0("(", paste(deparse(b, width.cutoff = 500L), collapse = " "), ")"), ""),
           vapply(specials, deparse1l, ""))
  f_new <- stats::as.formula(paste(paste(deparse(formula[[2]], width.cutoff = 500L), collapse = " "), "~",
                                   paste(rhs, collapse = " + ")), env = environment(formula))
  # the variables of the smooths (covariates, factors, 'by' variables): never
  # centred or scaled (a numeric 'by' multiplies the smooth)
  vars <- c(vars, unlist(lapply(rcurves, `[[`, "vars")))
  list(formula = f_new, groups = groups, fixed = unique(fixed), derived = derived, vars = unique(vars),
       rcurves = rcurves)
}

## The derived 'by' factors (cells of several factors) for data or new data.
pm_sm_derive <- function(data, derived, levels = NULL) {
  for (nm in names(derived)) {
    lab <- do.call(paste, c(lapply(data[derived[[nm]]], as.character), sep = ":"))
    data[[nm]] <- factor(lab, levels = if (is.null(levels[[nm]])) sort(unique(lab)) else levels[[nm]])
  }
  data
}

## Set up the smooths on 'data' and split them into fixed and random parts.
## 'sp': pm_sm_parse(). Returns list(blocks (named list of sparse n x q
## matrices, one variance each), fixed (character: fixed terms to add),
## smooths (what prediction needs: the mgcv smooth objects without their
## data, for each block the smooths it comes from with the linear maps, the
## derived 'by' factors), vars (variables of the smooths)).
pm_sm_build <- function(sp, data) {
  gd <- pm_sm_derive(data, sp$derived)
  lev <- lapply(sp$derived, function(x) NULL)
  for (nm in names(sp$derived)) lev[[nm]] <- levels(gd[[nm]])
  gd$.pm_y <- 0
  sms <- list()
  for (g in sp$groups) {
    # in the formula's environment: objects named in s() (xt = list(knots = kn)) are found there
    G <- mgcv::gam(stats::as.formula(paste(".pm_y ~", paste(g, collapse = " + ")), env = environment(sp$formula)),
                   data = gd, fit = FALSE)
    sms <- c(sms, G$smooth)
  }
  blocks <- list()
  maps <- list()
  fixed <- sp$fixed
  id_of <- character(0)
  for (i in seq_along(sms)) {
    sm <- sms[[i]]
    Xd <- mgcv::PredictMat(sm, gd)
    sm$X <- Xd
    r <- mgcv::smooth2random(sm, names(gd), type = 2)
    base <- pm_sm_name(sm$label, names(data))
    nb <- length(r$rand)
    for (j in seq_len(nb)) {
      R <- as.matrix(r$rand[[j]])
      nm <- if (nb == 1) base else paste0(base, "_", j)
      blocks[[nm]] <- Matrix::Matrix(R, sparse = TRUE)
      maps[[nm]] <- list(list(sm = i, A = pm_sm_map(sm, Xd, R, r)))
      id_of[nm] <- if (!is.null(sm$id) && nb == 1) as.character(sm$id) else NA_character_
    }
    if (NCOL(r$Xf)) {
      if (!is.null(sm$id) && startsWith(as.character(sm$id), "pm:")) {
        # pieces of an expanded smooth: the fixed terms of the expansion must span them
        M <- stats::model.matrix(stats::as.formula(paste("~", paste(sp$fixed, collapse = " + "))), data)
        if (max(abs(qr.resid(qr(M), as.matrix(r$Xf)))) > 1e-8 * max(1, max(abs(r$Xf)))) {
          stop("The fixed part of the smooth ", sm$label, " is not spanned by its fixed terms.", call. = FALSE)
        }
      } else {
        fixed <- union(fixed, pm_sm_fixed_terms(sm, as.matrix(r$Xf), data))
      }
    }
  }
  # smooths with the same id share one variance: one block (the deviation
  # curves of an expanded smooth: named by it)
  for (id in unique(stats::na.omit(id_of))) {
    nms <- names(id_of)[!is.na(id_of) & id_of == id]
    own <- startsWith(id, "pm:")
    if (length(nms) < 2 && !own) next
    new_nm <- if (own) sub("^pm:", "", id) else paste0(nms[1], "_id")
    blk <- if (length(nms) > 1) do.call(cbind, blocks[nms]) else blocks[[nms]]
    mp <- unlist(maps[nms], recursive = FALSE)
    blocks[nms] <- NULL
    maps[nms] <- NULL
    blocks[[new_nm]] <- blk
    maps[[new_nm]] <- mp
  }
  # random curves (s(x) | g): the wiggly part of the smooth, a copy per level
  # of g, one variance; random effects, so not part of the population curves
  rcurves <- list()
  for (rc in sp$rcurves) {
    G <- mgcv::gam(stats::as.formula(paste(".pm_y ~", rc$call), env = environment(sp$formula)), data = gd, fit = FALSE)
    sm <- G$smooth[[1]]
    Xd <- mgcv::PredictMat(sm, gd)
    sm$X <- Xd
    r <- mgcv::smooth2random(sm, names(gd), type = 2)
    if (length(r$rand) != 1) stop("A random curve (s(x) | g) needs a smooth with one penalty (e.g. bs = \"os\"): ",
                                  rc$call, call. = FALSE)
    R <- as.matrix(r$rand[[1]])
    g <- droplevels(factor(eval(str2lang(rc$group), gd)))
    n <- nrow(R)
    q <- ncol(R)
    blk <- Matrix::drop0(Matrix::sparseMatrix(i = rep(seq_len(n), q), j = (as.integer(g) - 1L)*q + rep(seq_len(q), each = n),
                                              x = as.vector(R), dims = c(n, q*nlevels(g))))
    nm <- pm_sm_name(paste0(sm$label, "|", rc$group), c(names(data), names(blocks)))
    blocks[[nm]] <- blk
    sm$X <- NULL
    rcurves[[nm]] <- list(sm = sm, A = pm_sm_map(sm, Xd, R, r), group = rc$group, levels = levels(g))
  }
  # each block scaled as a random intercept (sum of squares n): the variance
  # is free, so the model is the same, but lme4's start (theta = 1: the block
  # adds a variance like the residual's) is then near the optimum whatever the
  # covariate's units (mgcv's own scaling left some fits at a worse local
  # optimum)
  n_obs <- nrow(data)
  for (nm in names(blocks)) {
    sc <- sqrt(sum(blocks[[nm]]^2)/n_obs)
    if (!is.finite(sc) || sc <= 0) next
    blocks[[nm]] <- blocks[[nm]]/sc
    if (!is.null(maps[[nm]])) maps[[nm]] <- lapply(maps[[nm]], function(m) { m$A <- m$A/sc; m })
    if (!is.null(rcurves[[nm]])) rcurves[[nm]]$A <- rcurves[[nm]]$A/sc
  }
  for (k in seq_along(sms)) sms[[k]]$X <- NULL
  vars <- unique(c(setdiff(unlist(lapply(sms, function(sm) c(sm$term, if (sm$by != "NA") sm$by))),
                           names(sp$derived)), unlist(sp$derived), unlist(lapply(sp$rcurves, `[[`, "vars"))))
  list(blocks = blocks, fixed = unique(fixed), vars = vars,
       smooths = structure(list(sms = sms, maps = maps, labels = vapply(sms, `[[`, "", "label"),
                                derived = sp$derived, levels = lev, rcurves = rcurves), class = "pm_smooths"))
}

## A syntactic name for a smooth block, not a column of the data.
pm_sm_name <- function(label, taken) {
  nm <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", label))
  while (nm %in% taken) nm <- paste0(nm, "_sm")
  nm
}

## The fixed terms that span the null space of a smooth (its fixed columns
## Xf on the data, from smooth2random()), with the intercept: the variables
## of the smooth (x1 + x2), then their interaction (x1*x2), then squares;
## times the 'by' variable (f + f:x for a factor, z + z:x for a numeric one).
pm_sm_fixed_terms <- function(sm, Xf, data) {
  v <- sm$term
  by <- if (sm$by != "NA") sm$by
  cands <- list(v)
  if (length(v) > 1) cands <- c(cands, list(c(v, attr(stats::terms(stats::as.formula(paste("~", paste(v, collapse = "*")))),
                                                     "term.labels"))))
  cands <- c(cands, list(c(cands[[length(cands)]], paste0("I(", v, "^2)"))))
  spans <- function(trm) {
    f <- stats::as.formula(paste("~", paste(trm, collapse = " + ")))
    M <- stats::model.matrix(f, data)
    max(abs(qr.resid(qr(M), Xf))) <= 1e-8 * max(1, max(abs(Xf)))
  }
  for (cand in cands) {
    trm <- if (is.null(by)) cand else c(by, paste0(by, ":", cand))
    if (spans(trm)) return(trm)
  }
  stop("Can't find the fixed part of the smooth ", sm$label, "; write it in the formula.", call. = FALSE)
}

## The random-effect design of the smooth blocks for new data (a named list
## of sparse matrices, in the order of the blocks).
pm_sm_predict <- function(smooths, newdata, blocks = names(smooths$maps)) {
  if (length(smooths$derived)) newdata <- pm_sm_derive(newdata, smooths$derived, smooths$levels)
  P <- list()
  out <- lapply(blocks, function(b) {
    parts <- lapply(smooths$maps[[b]], function(m) {
      k <- as.character(m$sm)
      if (is.null(P[[k]])) P[[k]] <<- mgcv::PredictMat(smooths$sms[[m$sm]], newdata)
      P[[k]] %*% m$A
    })
    Matrix::Matrix(do.call(cbind, parts), sparse = TRUE)
  })
  stats::setNames(out, blocks)
}

## The linear map A from a smooth's design Xd (PredictMat() on the data) to
## one of its random parts R (Xd A = R), used again for prediction. For a
## smooth with one penalty it is smooth2random()'s own transformation (the
## eigenvectors of the penalty times 1/sqrt(eigenvalues)), exact even when
## Xd is not of full rank on the data (fewer distinct covariate values than
## basis columns). Otherwise it is found by regression; a factor smooth (bs =
## "fs") has a copy of its base design per level, and when a level has too
## few distinct values the map is built from the base design (the same for
## every level).
pm_sm_map <- function(sm, Xd, R, r = NULL) {
  ok <- function(A) max(abs(as.matrix(Xd %*% A) - R)) <= 1e-8*max(1, max(abs(R)))
  if (!is.null(r$trans.U) && length(r$rand) == 1 && ncol(r$trans.U) == ncol(Xd)) {
    A <- t(t(r$trans.U)*r$trans.D)[, seq_len(ncol(R)), drop = FALSE]
    if (ok(A)) return(A)
  }
  A <- tryCatch(qr.solve(Xd, R), error = function(e) NULL)
  if (!is.null(A)) return(A)
  if (inherits(sm, "fs.interaction") && !is.null(sm$Xb)) {
    nl <- length(sm$flev)
    q <- ncol(R) %/% nl
    Rb <- Reduce(`+`, lapply(seq_len(nl), function(j) R[, (j - 1)*q + seq_len(q), drop = FALSE]))
    B <- tryCatch(qr.solve(sm$Xb, Rb), error = function(e) NULL)
    if (!is.null(B)) {
      A <- kronecker(diag(nl), B)
      if (ok(A)) return(A)
    }
  }
  stop("The smooth ", sm$label, " cannot be set up: its basis is not of full rank on the data ",
       "(too few distinct covariate values for 'k'?).", call. = FALSE)
}

## ---- The own bases as mgcv smooths ----
##
## bs = "os": the O'Sullivan penalised spline of Wand and Ormerod (2008),
## as smZ(x, k) and Harezlak, Ruppert and Wand (2018): k interior knots at
## quantiles of the unique x (default 6, as smZ()), the range widened by
## 1% (xt = list(knots = <interior knots>, range = c(lo, hi)) sets them).
## bs = "ztps": the thin-plate basis of Ztps() for two covariates (k knots
## per axis; xt = list(knots = <two-column matrix>) sets the knots). Both give the design [1, x, Z] with the identity penalty
## on Z and no constraint (the constant and x are the fixed part), so the
## random effects have exactly smZ()'s design. (Registered as S3 methods of
## mgcv's generics in NAMESPACE.)

smooth.construct.os.smooth.spec <- function(object, data, knots) {
  x <- data[[object$term]]
  k <- if (object$bs.dim < 0) 6 else object$bs.dim
  k <- min(k, length(unique(x)) - 2)
  xt <- if (is.list(object$xt)) object$xt else list()
  intKnots <- if (!is.null(xt$knots)) sort(as.numeric(xt$knots)) else if (!is.null(knots[[object$term]]))
    knots[[object$term]] else stats::quantile(unique(x), seq(0, 1, length = k + 2)[-c(1, k + 2)])
  range.x <- if (!is.null(xt$range)) as.numeric(xt$range) else c(1.01*min(x) - 0.01*max(x), 1.01*max(x) - 0.01*min(x))
  Z <- ZOSull(x, intKnots = intKnots, range.x = range.x)
  q <- ncol(Z)
  object$X <- cbind(1, x, Z)
  object$S <- list(diag(c(0, 0, rep(1, q))))
  object$rank <- q
  object$null.space.dim <- 2
  object$C <- matrix(0, 0, ncol(object$X))
  object$bs.dim <- ncol(object$X)
  object$df <- ncol(object$X)
  object$intKnots <- intKnots
  object$range.x <- range.x
  class(object) <- c("os.smooth", "mgcv.smooth")
  object
}

Predict.matrix.os.smooth <- function(object, data) {
  x <- data[[object$term]]
  cbind(1, x, ZOSull(x, intKnots = object$intKnots, range.x = object$range.x))
}

smooth.construct.ztps.smooth.spec <- function(object, data, knots) {
  if (length(object$term) != 2) stop("bs = \"ztps\" is for two covariates.", call. = FALSE)
  x <- as.matrix(as.data.frame(data[object$term]))
  k <- if (object$bs.dim < 0) 10 else object$bs.dim
  kn0 <- if (is.list(object$xt) && !is.null(object$xt$knots)) as.matrix(object$xt$knots)
  Z <- if (is.null(kn0)) Ztps(x, k = k) else Ztps(x, knots = kn0)
  kn <- attr(Z, "knots")
  Z <- unclass(as.matrix(Z))
  q <- ncol(Z)
  object$X <- cbind(1, x, Z)
  object$S <- list(diag(c(0, 0, 0, rep(1, q))))
  object$rank <- q
  object$null.space.dim <- 3
  object$C <- matrix(0, 0, ncol(object$X))
  object$bs.dim <- ncol(object$X)
  object$df <- ncol(object$X)
  object$ztps_k <- k
  object$ztps_knots <- kn
  class(object) <- c("ztps.smooth", "mgcv.smooth")
  object
}

Predict.matrix.ztps.smooth <- function(object, data) {
  x <- as.matrix(as.data.frame(data[object$term]))
  Z <- unclass(as.matrix(Ztps(x, k = object$ztps_k, knots = object$ztps_knots)))
  cbind(1, x, Z)
}
