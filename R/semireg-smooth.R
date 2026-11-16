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
## random terms and without the smooths), smooths (character: the smooth
## calls), vars (their variables)) or NULL when there is none.
pm_sm_parse <- function(formula) {
  bars <- pm_findbars(formula)
  f_fixed <- pm_nobars(formula)
  is_sm <- function(e) is.call(e) && as.character(e[[1]])[1] %in% c("s", "te", "t2", "ti")
  rhs_terms <- attr(stats::terms(f_fixed, specials = c("s", "te", "t2", "ti")), "variables")[-1]
  sm_calls <- Filter(is_sm, as.list(rhs_terms))
  if (!length(sm_calls)) return(NULL)
  sm_txt <- vapply(sm_calls, function(e) {
    # tensor products: mgcv's pre-fit design of t2() is not in the span of its
    # prediction matrix (mgcv 1.9-1), so its split into fixed and random parts
    # cannot be checked yet (dev/review-semireg.md); te() and ti() cannot be
    # written as independent random effects at all
    if (as.character(e[[1]]) %in% c("te", "t2", "ti")) {
      stop("Tensor-product smooths (te(), t2(), ti()) are not available in semireg() yet; use s(x1, x2) ",
           "(isotropic thin plate) or s(x1, x2, bs = \"ztps\").", call. = FALSE)
    }
    paste(deparse(e, width.cutoff = 500L), collapse = " ")
  }, "")
  # a smooth may only enter as a main term (not s(x):A)
  labs <- attr(stats::terms(f_fixed), "term.labels")
  sm_lab <- vapply(sm_calls, function(e) paste(deparse(e, width.cutoff = 500L), collapse = " "), "")
  bad <- labs[!labs %in% sm_lab & vapply(labs, function(l) any(vapply(sm_lab, grepl, logical(1), x = l,
                                                                         fixed = TRUE)), logical(1))]
  if (length(bad)) stop("A smooth can only be a term of its own: ", paste(bad, collapse = ", "), call. = FALSE)
  ig <- mgcv::interpret.gam(f_fixed)
  pf <- ig$pf
  # the variables of the smooths (covariates and 'by' variables): never
  # centred or scaled (a numeric 'by' multiplies the smooth)
  vars <- unique(unlist(lapply(ig$smooth.spec, function(sp) c(sp$term, if (sp$by != "NA") sp$by))))
  rhs <- c(paste(deparse(pf[[3]], width.cutoff = 500L), collapse = " "),
           vapply(bars, function(b) paste0("(", paste(deparse(b, width.cutoff = 500L), collapse = " "), ")"), ""))
  f_new <- stats::as.formula(paste(paste(deparse(formula[[2]], width.cutoff = 500L), collapse = " "), "~",
                                   paste(rhs, collapse = " + ")), env = environment(formula))
  list(formula = f_new, smooths = sm_txt, vars = vars)
}

## Set up the smooths on 'data' and split them into fixed and random parts.
## Returns list(blocks (named list of sparse n x q matrices, one variance
## each), fixed (character: fixed terms to add), smooths (what prediction
## needs: the mgcv smooth objects without their data, and for each block
## the smooths it comes from with the linear maps), vars (variables of the
## smooths)).
pm_sm_build <- function(smooths, data) {
  gd <- data
  gd$.pm_y <- 0
  fg <- stats::as.formula(paste(".pm_y ~", paste(smooths, collapse = " + ")))
  G <- mgcv::gam(fg, data = gd, fit = FALSE)
  sms <- G$smooth
  vars <- unique(unlist(lapply(sms, function(sm) c(sm$term, if (sm$by != "NA") sm$by))))
  blocks <- list()
  maps <- list()
  fixed <- character(0)
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
      maps[[nm]] <- list(list(sm = i, A = qr.solve(Xd, R)))
      id_of[nm] <- if (!is.null(sm$id) && nb == 1) as.character(sm$id) else NA_character_
    }
    if (NCOL(r$Xf)) fixed <- union(fixed, pm_sm_fixed_terms(sm, as.matrix(r$Xf), data))
  }
  # smooths with the same id share one variance: one block
  for (id in unique(stats::na.omit(id_of))) {
    nms <- names(id_of)[!is.na(id_of) & id_of == id]
    if (length(nms) < 2) next
    new_nm <- paste0(nms[1], "_id")
    blocks[[new_nm]] <- do.call(cbind, blocks[nms])
    maps[[new_nm]] <- unlist(maps[nms], recursive = FALSE)
    blocks[nms] <- NULL
    maps[nms] <- NULL
  }
  for (k in seq_along(sms)) sms[[k]]$X <- NULL
  list(blocks = blocks, fixed = fixed, vars = vars,
       smooths = structure(list(sms = sms, maps = maps, labels = vapply(sms, `[[`, "", "label")),
                           class = "pm_smooths"))
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

## ---- The own bases as mgcv smooths ----
##
## bs = "os": the O'Sullivan penalised spline of Wand and Ormerod (2008),
## as smZ(x, k) and Harezlak, Ruppert and Wand (2018): k interior knots at
## quantiles of the unique x (default 6, as smZ()), the range widened by
## 1%. bs = "ztps": the thin-plate basis of Ztps() for two covariates (k
## knots per axis). Both give the design [1, x, Z] with the identity penalty
## on Z and no constraint (the constant and x are the fixed part), so the
## random effects have exactly smZ()'s design. (Registered as S3 methods of
## mgcv's generics in NAMESPACE.)

smooth.construct.os.smooth.spec <- function(object, data, knots) {
  x <- data[[object$term]]
  k <- if (object$bs.dim < 0) 6 else object$bs.dim
  k <- min(k, length(unique(x)) - 2)
  intKnots <- if (!is.null(knots[[object$term]])) knots[[object$term]] else
    stats::quantile(unique(x), seq(0, 1, length = k + 2)[-c(1, k + 2)])
  range.x <- c(1.01*min(x) - 0.01*max(x), 1.01*max(x) - 0.01*min(x))
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
  Z <- Ztps(x, k = k)
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
