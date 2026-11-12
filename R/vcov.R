## The 'vcov' argument of predictmeans(), contrastmeans(), pm_slopes() and
## pm_joint(): a covariance of the coefficients other than the model's.
##
##   NULL                  the model's (with vcov_adjust for lmer), as before;
##   "HC0" ... "HC3"       heteroscedasticity-consistent (sandwich) covariance
##                         for lm, aov and glm, equal to sandwich::vcovHC();
##   ~ cluster             cluster-robust (CR1) covariance for lm, aov and glm,
##                         equal to sandwich::vcovCL() (type "HC1": the factor
##                         G/(G-1) (n-1)/(n-p));
##   a matrix              used as given (rows and columns named by the
##                         coefficients, or in their order);
##   a function            called with the model, e.g. a function from
##                         clubSandwich for mixed models.

## Sandwich covariance of an lm or glm fit: (X'WX)^-1 M (X'WX)^-1, M the sum
## of outer products of the estimating functions w u x (u working residuals,
## w working weights), each scaled for HC1-HC3, or summed within clusters.
pm_vcov_sandwich <- function(model, type = "HC3", cluster = NULL) {
  if (!inherits(model, "lm")) {
    stop("vcov = \"", type, "\" and vcov = ~cluster are for lm, aov and glm models; for other models give ",
         "vcov a covariance matrix or a function (for mixed models e.g. from the clubSandwich package).",
         call. = FALSE)
  }
  b <- stats::coef(model)
  ok <- !is.na(b)
  X <- stats::model.matrix(model)[, ok, drop = FALSE]
  if (inherits(model, "glm")) {
    w <- model$weights
    u <- stats::residuals(model, type = "working")
  } else {
    w <- stats::weights(model)
    if (is.null(w)) w <- rep(1, nrow(X))
    u <- stats::residuals(model)
  }
  use <- w > 0
  X <- X[use, , drop = FALSE]
  w <- w[use]
  u <- u[use]
  n <- nrow(X)
  p <- ncol(X)
  Xw <- X * sqrt(w)
  B <- solve(crossprod(Xw))
  ef <- X * (w * u)                                # estimating functions, one row per observation
  if (!is.null(cluster)) {
    cluster <- cluster[use]
    G <- length(unique(cluster))
    if (G < 2) stop("vcov = ~cluster needs at least 2 clusters.", call. = FALSE)
    ef <- rowsum(ef, cluster, reorder = FALSE)
    adj <- G / (G - 1) * (n - 1) / (n - p)
  } else {
    h <- rowSums((Xw %*% B) * Xw)
    adj <- 1
    if (type == "HC1") adj <- n / (n - p)
    if (type == "HC2") ef <- ef / sqrt(1 - h)
    if (type == "HC3") ef <- ef / (1 - h)
  }
  V <- adj * B %*% crossprod(ef) %*% B
  dimnames(V) <- list(names(b)[ok], names(b)[ok])
  V
}

## The cluster variable of vcov = ~ g, for the observations used in the fit.
pm_vcov_cluster <- function(model, f) {
  vars <- all.vars(f)
  if (length(vars) != 1) stop("vcov = ~cluster takes one variable, e.g. vcov = ~ block.", call. = FALSE)
  dat <- pm_get_data(model)
  if (!vars %in% names(dat)) stop("cluster variable '", vars, "' is not in the model's data.", call. = FALSE)
  rn <- rownames(stats::model.frame(model))
  cl <- if (!is.null(rownames(dat)) && all(rn %in% rownames(dat))) dat[rn, vars] else dat[[vars]]
  if (length(cl) != length(rn)) stop("cannot match the cluster variable to the model's observations.", call. = FALSE)
  if (anyNA(cl)) stop("the cluster variable has missing values.", call. = FALSE)
  cl
}

## Resolve 'vcov' to a matrix matching the model's estimable coefficients
## (names(mp$coef)); attribute "label" describes it for the settings table.
pm_vcov_resolve <- function(model, vcov, coef_names) {
  if (is.null(vcov)) return(NULL)
  if (is.logical(vcov)) {
    stop("'vcov' takes \"HC0\" to \"HC3\", ~cluster, a matrix or a function; for Kenward-Roger standard errors ",
         "use 'vcov_adjust'.", call. = FALSE)
  }
  if (inherits(vcov, "formula")) {
    V <- pm_vcov_sandwich(model, cluster = pm_vcov_cluster(model, vcov))
    label <- paste0("CR1, clusters: ", all.vars(vcov))
  } else if (is.character(vcov)) {
    type <- toupper(vcov[1])
    if (!type %in% c("HC0", "HC1", "HC2", "HC3")) {
      stop("vcov = \"", vcov[1], "\" is not available; use \"HC0\", \"HC1\", \"HC2\" or \"HC3\".", call. = FALSE)
    }
    V <- pm_vcov_sandwich(model, type)
    label <- type
  } else if (is.function(vcov)) {
    V <- vcov(model)
    label <- "user function"
  } else {
    V <- vcov
    label <- "user matrix"
  }
  V <- as.matrix(V)
  k <- length(coef_names)
  if (!is.null(rownames(V)) && !is.null(colnames(V)) && all(coef_names %in% rownames(V)) &&
      all(coef_names %in% colnames(V))) {
    V <- V[coef_names, coef_names, drop = FALSE]
  } else if (nrow(V) != k || ncol(V) != k) {
    stop("'vcov' must be a ", k, " x ", k, " matrix for the model's coefficients (",
         paste(utils::head(coef_names, 4), collapse = ", "), if (k > 4) ", ...", ").", call. = FALSE)
  } else {
    dimnames(V) <- list(coef_names, coef_names)
  }
  if (anyNA(V)) stop("'vcov' has missing values.", call. = FALSE)
  attr(V, "label") <- label
  V
}

## With a permlist the permutation statistics use their own covariance, so a
## 'vcov' is not used (for robust permutation tests: permmodels(robust = TRUE)).
pm_vcov_with_perm <- function(vcov, permlist) {
  if (is.null(vcov) || !pm_perm_given(permlist)) return(vcov)
  warning("'vcov' is not used with 'permlist'; for heteroscedasticity-robust permutation tests use ",
          "permmodels(..., robust = TRUE).", call. = FALSE)
  NULL
}
