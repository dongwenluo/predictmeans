#' Heritability on a Line-Mean (Entry-Mean) Basis
#'
#' The share of the variance of genotype means that is genetic, from a mixed model with random genotype
#' effects: the standard form for balanced trials and the mixed-model forms of Cullis et al. (2006), Piepho and
#' Möhring (2007) and Oakey et al. (2006) for unbalanced trials, incomplete blocks or correlated errors.
#'
#' @details
#' With \eqn{\sigma^2_g} the genotypic variance:
#' \describe{
#'   \item{\code{"standard"}}{\eqn{\sigma^2_g / (\sigma^2_g + \sum_t \sigma^2_t / n_t + \sigma^2_e / n)}, over
#'   the random terms \eqn{t} that contain the genotype (genotype x environment, ...), with \eqn{n_t} the number of
#'   their levels per genotype and \eqn{n} the number of plots per genotype (harmonic means over genotypes; give
#'   \code{n_env}, \code{n_rep} to set them). With one environment and \eqn{r} replicates,
#'   \eqn{\sigma^2_g/(\sigma^2_g + \sigma^2_e/r)}; with \eqn{e} environments,
#'   \eqn{\sigma^2_g/(\sigma^2_g + \sigma^2_{ge}/e + \sigma^2_e/(re))}. Exact for balanced data only.}
#'   \item{\code{"cullis"}}{\eqn{1 - \bar v_{BLUP}/(2\sigma^2_g)}, \eqn{\bar v_{BLUP}} the mean variance of a
#'   difference of two genotype BLUPs (prediction error variance, from the inverse of the mixed model
#'   equations, as \code{se_ranef(type = "prediction")}).}
#'   \item{\code{"piepho"}}{\eqn{\sigma^2_g / (\sigma^2_g + \bar v_{BLUE}/2)}, \eqn{\bar v_{BLUE}} the mean
#'   variance of a difference of two genotype BLUEs from the same model with genotype as a fixed effect (other
#'   random terms kept). By default (\code{piepho = "fixed"}) its variance components are those of the
#'   random-genotype fit, as Schmidt et al. (2019) do, so both terms of the ratio rest on the same estimates (no
#'   refit); \code{piepho = "refit"} estimates them again, as Piepho and Möhring (2007). On the oats alpha
#'   design of John and Williams (1995): 0.808 (Schmidt et al. 2019, Example 2) and 0.803 (Piepho and Möhring
#'   2007, Example 1). (\code{inti::H2cal(emmeans = TRUE)} and \pkg{heritable} take emmeans's Kenward-Roger
#'   adjusted covariance: 0.797 there.)}
#'   \item{\code{"oakey"}}{the mean of the non-zero eigenvalues of \eqn{I - C_{gg}/\sigma^2_g}, \eqn{C_{gg}}
#'   the prediction error covariance of the genotype BLUPs (the zero ones are directions confounded with fixed
#'   effects, as the mean of all genotypes with the intercept). With independent genotypes it equals
#'   \code{"cullis"} when the intercept is the only fixed effect the genotypes are confounded with; with genotypes
#'   nested in a fixed factor (families in provenances) the two differ (a note says so), and with a relationship
#'   matrix.}
#'   \item{\code{"reliability"}}{the mean reliability of the genotype BLUPs,
#'   \eqn{\bar r^2 = m^{-1}\sum_i (1 - C_{ii}/\sigma^2_g)} (Schmidt et al. 2019, Eq. 44); it ignores the
#'   covariances of the prediction errors, so it is a little below the others for small \eqn{m}.}
#'   \item{\code{"simulated"}}{the expected squared correlation between true genotypic values and their BLUPs,
#'   simulated from their joint covariance (Piepho and Möhring 2007; Schmidt et al. 2019, Eq. 41), with
#'   \code{nsim} draws.}
#'   \item{\code{"family"}}{narrow-sense heritability of a single tree (plant, animal) from families:
#'   \eqn{(\sigma^2_f/\rho) / (\sigma^2_f + \sum_t \sigma^2_t + \sigma^2_e)}, \eqn{\sigma^2_f} the family variance
#'   (\code{genotype} is the family factor), \eqn{\rho} the \code{relationship} among members of a family (1/4
#'   for half-sibs, 1/2 for full-sibs; often taken as 1/3 for open-pollinated progeny, which are partly full-sibs
#'   or selfs), and \eqn{\sigma^2_t} the random terms with the family in them (family x block, i.e. plot; family x
#'   site). The other design terms (blocks, rows, columns) are left out of the phenotypic variance. Can exceed
#'   1 when \eqn{\rho} is too small for the material. Family-mean heritability is the line-mean forms above
#'   with the family factor as \code{genotype}.}
#'   \item{\code{"repeatability"}}{\eqn{\sigma^2_s / (\sigma^2_s + \sum_t \sigma^2_t + \sigma^2_e)} with
#'   \code{genotype} the subject that has repeated records: the share of the variance of a single record due to
#'   the subject.}
#' }
#' The four line-mean forms agree for balanced data with independent genotypes; with imbalance, incomplete blocks or
#' spatial errors they differ, and Schmidt et al. (2019) recommend the mixed-model forms. Narrow-sense
#' heritability from a pedigree or genomic relationship matrix needs a model with that matrix, which lme4
#' cannot fit: \pkg{heritable} computes these forms (and the per-genotype \eqn{H^2_\Delta} of Schmidt et
#' al. 2019) from asreml fits with a relationship matrix. Its standard, Cullis and Oakey values for lmer fits
#' are those of this function (tested); this function adds standard errors, the family and repeatability
#' forms, and advice.
#' @param model An \code{lmer} fit, or a gaussian \code{semireg} fit, with a random genotype term
#'   \code{(1 | genotype)}. A genotypic variance estimated as 0 gives 0 for every form, with a note.
#' @param genotype Name of the genotype (line, entry, family) factor.
#' @param method \code{"cullis"} (default), \code{"piepho"}, \code{"standard"}, \code{"oakey"},
#'   \code{"reliability"}, \code{"simulated"}, \code{"family"}, \code{"repeatability"}, or \code{"all"} (the
#'   line-mean forms: standard to simulated).
#' @param piepho For \code{"piepho"}: \code{"fixed"} (default) or \code{"refit"}, see Details.
#' @param nsim For \code{"simulated"}: the number of simulated draws.
#' @param n_env,n_rep For \code{"standard"}, a model with one genotype x environment term: the number of
#'   environments and of replicates per environment (\code{n_rep} alone: replicates, with the environments from
#'   the data); by default both from the data (harmonic means), the only way with several such terms.
#' @param relationship For \code{"family"}: the coefficient of relationship among members of a family.
#' @param within_variance,within_n For \code{"family"} and \code{"repeatability"} from plot means: the variance
#'   among plants within plots and the number of plants per plot (numbers, or names of columns of the data:
#'   the mean of the variances, the harmonic mean of the numbers). The single-plant variance is then
#'   \eqn{\sigma^2_f + \sum_t \sigma^2_t + (\sigma^2_e - \sigma^2_w/n) + \sigma^2_w}, as Williams et al. (2023,
#'   Example 6.2), with \eqn{\sigma^2_e - \sigma^2_w/n} the plot variance; \eqn{\sigma^2_w} is taken as known in
#'   the standard error.
#' @param nboot For the mixed-model forms (Cullis to simulated): parametric bootstrap draws for their
#'   standard error and percentile limits (each draw refits the model; 0, the default, gives none). The other
#'   forms get theirs by the delta method from the covariance of the variance parameters, as \code{varcomp()}.
#' @param ncore Cores for the bootstrap; results are the same for any number.
#' @param seed Seed for the simulated and bootstrap draws.
#' @param conf_level Confidence level of the limits (logit-normal for the delta method, so they stay in (0, 1)).
#' @return A data frame of class \code{"pm_heritability"}: \code{method}, \code{estimate}, \code{std.error},
#'   \code{conf.low}, \code{conf.high}; \code{attr(, "components")} holds the variance components and the mean
#'   variances of differences used.
#' @references
#' Cullis, B. R., Smith, A. B. and Coombes, N. E. (2006). On the design of early generation variety trials
#' with correlated data. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, 11, 381-393.
#'
#' Oakey, H., Verbyla, A., Pitchford, W., Cullis, B. and Kuchel, H. (2006). Joint modeling of additive and
#' non-additive genetic line effects in single field trials. \emph{Theoretical and Applied Genetics}, 113,
#' 809-819.
#'
#' Piepho, H.-P. and Möhring, J. (2007). Computing heritability and selection response from unbalanced plant
#' breeding trials. \emph{Genetics}, 177, 1881-1888.
#'
#' Schmidt, P., Hartung, J., Bennewitz, J. and Piepho, H.-P. (2019). Heritability in plant breeding on a
#' genotype-difference basis. \emph{Genetics}, 212(4), 991-1008. \doi{10.1534/genetics.119.302134}
#' @seealso \code{\link{varcomp}}, \code{\link{se_ranef}}
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lmer(yield ~ nitro + (1 | Block) + (1 | Variety), data = Oats)
#' heritability(fm, "Variety", method = "all")
#' \donttest{
#' heritability(fm, "Variety", nboot = 199, seed = 1)   # SE of the Cullis form by bootstrap
#' }
#' @export
heritability <- function(model, genotype,
                         method = c("cullis", "piepho", "standard", "oakey", "reliability", "simulated", "family",
                                    "repeatability", "all"),
                         relationship = NULL, within_variance = NULL, within_n = NULL,
                         n_env = NULL, n_rep = NULL, piepho = c("fixed", "refit"),
                         nsim = 10000, nboot = 0, ncore = 1, seed = NULL, conf_level = 0.95) {
  method <- match.arg(method)
  piepho <- match.arg(piepho)
  methods <- if (method == "all") c("standard", "cullis", "piepho", "oakey", "reliability", "simulated") else method
  if (!is.null(seed)) set.seed(seed)
  if (method == "family" && !(is.numeric(relationship) && length(relationship) == 1 && relationship > 0 &&
                              relationship <= 1)) {
    stop("method = \"family\" needs 'relationship', the coefficient of relationship among members of a family ",
         "(1/4 half-sibs, 1/2 full-sibs; often 1/3 for open-pollinated progeny).", call. = FALSE)
  }
  alpha <- pm_conf_alpha(conf_level)
  user_vars <- tryCatch(names(pm_get_data(model)), error = function(e) NULL)   # (not semireg's smooth groups)
  within <- pm_h2_within(model, within_variance, within_n)
  model <- pm_fit(model)
  if (!inherits(model, "lmerMod")) {
    stop("heritability() handles lmer and gaussian semireg fits (random genotype effects); a ", class(model)[1],
         " model is not yet covered.", call. = FALSE)
  }
  need_C <- any(methods %in% c("cullis", "oakey", "reliability", "simulated"))
  comp <- pm_h2_components(model, genotype, n_env, n_rep, need_C = need_C, user_vars = user_vars)
  comp$piepho <- piepho
  comp$nsim <- nsim
  comp$within <- within
  est <- pm_h2_point(model, comp, setdiff(methods, "simulated"), relationship)
  mc_se <- NULL
  if ("simulated" %in% methods) {
    sim <- if (comp$s2g > 0) pm_h2_sim(comp$C, comp$s2g, nsim) else structure(0, mc_se = 0)
    est <- c(est, simulated = as.numeric(sim))[methods]
    mc_se <- attr(sim, "mc_se")
  }
  # SE and limits: the delta method on the variance parameters for the forms
  # that are functions of the components; a parametric bootstrap for the
  # forms that also depend on the prediction errors (nboot > 0)
  vfun <- pm_h2_vfun(comp, relationship)
  boot_m <- intersect(methods, c("cullis", "piepho", "oakey", "reliability", "simulated"))
  bt <- if (length(boot_m) && nboot > 0) pm_h2_boot(model, comp, boot_m, nboot, ncore, seed)
  rows <- lapply(methods, function(m) {
    if (m %in% names(vfun)) {
      se <- if (comp$s2g > 0) pm_h2_delta(model, vfun[[m]]) else NA_real_
      # the family form ranges over (0, 1/relationship): limits on that range
      lim <- pm_h2_limits(est[[m]], se, alpha, scale = if (m == "family") relationship else 1)
    } else if (!is.null(bt)) {
      draws <- bt[, m]
      se <- stats::sd(draws, na.rm = TRUE)
      lim <- unname(stats::quantile(draws, c(alpha/2, 1 - alpha/2), na.rm = TRUE))
    } else {
      se <- NA_real_
      lim <- c(NA_real_, NA_real_)
    }
    data.frame(method = m, estimate = est[[m]], std.error = se, conf.low = lim[1], conf.high = lim[2])
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- structure(out, class = c("pm_heritability", "data.frame"), conf_level = conf_level,
                   components = comp[c("s2g", "terms", "s2e", "n_terms", "n_plots", "within")], relationship = relationship,
                   genotype = genotype, mc_se = mc_se, nsim = nsim,
                   nboot = if (!is.null(bt)) sum(stats::complete.cases(bt)) else 0)
  pm_add_advice(out, tryCatch(pm_heritability_advice(out, comp, nboot), error = function(e) NULL))
}

## The estimates of the requested forms from one fit (also each bootstrap refit).
pm_h2_point <- function(model, comp, methods, relationship = NULL, y = NULL) {
  # a genotypic variance estimated as 0 (a boundary fit): every form is 0
  if (comp$s2g <= 0) return(stats::setNames(rep(0, length(methods)), methods))
  vf <- pm_h2_vfun(comp, relationship)
  v <- pm_h2_vc(comp)
  vapply(methods, function(m) switch(m,
    cullis = 1 - pm_h2_vdiff(comp$C)/(2*comp$s2g),
    piepho = comp$s2g / (comp$s2g + (if (identical(comp$piepho, "refit")) pm_h2_vblue(model, comp$genotype, y) else
      pm_h2_vblue_fixed(model, comp$genotype))/2),
    reliability = mean(1 - diag(comp$C)/comp$s2g),
    simulated = as.numeric(pm_h2_sim(comp$C, comp$s2g, comp$nsim)),
    oakey = pm_h2_oakey(comp$C, comp$s2g),
    vf[[m]](v)), numeric(1))
}

## The variance components of the fit, named as pm_vc_mer_fun() names them.
pm_h2_vc <- function(comp) {
  stats::setNames(c(comp$s2g, comp$terms, comp$s2e),
                  c(paste0(comp$genotype, ".(Intercept)"), pm_h2_tnames(comp$with_g), "residual"))
}

## The forms that are functions of the variance components alone (v, named as
## pm_vc_mer_fun()): standard (line mean), family (single tree, narrow sense:
## the family variance is relationship x additive variance) and repeatability
## (single record). The terms with the genotype in them (genotype x
## environment, family x block = plot) are in the variance of a single record;
## the other design terms (blocks, rows, columns) are not.
pm_h2_vfun <- function(comp, relationship = NULL) {
  g <- paste0(comp$genotype, ".(Intercept)")
  tt <- pm_h2_tnames(comp$with_g)
  # plot means: the residual holds s2w / n of the plants within plots; a
  # single plant has all of s2w (s2w: the within-plot variance, n: plants per plot)
  w <- comp$within
  extra <- if (is.null(w)) 0 else w$s2w * (1 - 1/w$n)
  single <- function(v) v[[g]] + sum(v[tt]) + v[["residual"]] + extra
  list(standard = function(v) v[[g]] / (v[[g]] + sum(v[tt] / comp$n_terms) + v[["residual"]] / comp$n_plots),
       family = function(v) (v[[g]] / relationship) / single(v),
       repeatability = function(v) v[[g]] / single(v))
}

## The pieces of every form, from one fit: the genotypic variance, the other
## random terms that contain the genotype with their numbers of levels per
## genotype, the residual, the plots per genotype, and the prediction error
## covariance of the genotype BLUPs (pm_ranef_cov(), R/ranef-cov.R).
pm_h2_components <- function(model, genotype, n_env = NULL, n_rep = NULL, need_C = TRUE, user_vars = NULL) {
  if (stats::family(model)$family != "gaussian") stop("heritability() needs a gaussian lmer fit.", call. = FALSE)
  vc <- lme4::VarCorr(model)
  # the genotype's random intercept (lme4 names a second term of the same
  # factor, e.g. a random slope, "gen.1")
  gi <- which(sub("\\.[0-9]+$", "", names(vc)) == genotype &
                vapply(vc, function(x) identical(colnames(x), "(Intercept)"), logical(1)))
  if (!length(gi)) {
    stop("heritability() needs a random intercept for '", genotype, "', (1 | ", genotype, "), in the model.",
         call. = FALSE)
  }
  s2g <- max(vc[[gi[1]]][1, 1], 0)
  grp <- strsplit(names(vc), ":", fixed = TRUE)
  with_g <- names(vc)[vapply(grp, function(g) all(strsplit(genotype, ":", fixed = TRUE)[[1]] %in% g) &&
                               length(g) > length(strsplit(genotype, ":", fixed = TRUE)[[1]]), logical(1))]
  dat <- stats::model.frame(model)
  gcol <- pm_h2_gcol(dat, genotype)
  per_g <- function(vars) {             # levels of 'vars' per genotype (harmonic mean)
    k <- tapply(do.call(paste, c(dat[vars], sep = "\r")), gcol, function(x) length(unique(x)))
    k <- k[!is.na(k)]
    length(k) / sum(1/k)
  }
  if ((!is.null(n_env) || !is.null(n_rep)) && length(with_g) > 1) {
    stop("'n_env' and 'n_rep' are for a model with one genotype x environment term; this one has ",
         paste(with_g, collapse = ", "), ": leave them out to take the numbers from the data.", call. = FALSE)
  }
  if (!is.null(n_env) && !length(with_g)) {
    stop("'n_env' needs a genotype x environment term in the model.", call. = FALSE)
  }
  n_terms <- vapply(with_g, function(t) {
    other <- setdiff(strsplit(t, ":", fixed = TRUE)[[1]], strsplit(genotype, ":", fixed = TRUE)[[1]])
    if (!is.null(n_env)) n_env else per_g(other)
  }, numeric(1))
  kp <- table(gcol); kp <- kp[kp > 0]
  n_plots <- if (!is.null(n_rep)) n_rep * (if (length(n_terms)) n_terms[[1]] else 1) else length(kp) / sum(1/kp)
  s2t <- vapply(with_g, function(t) vc[[t]][1, 1], numeric(1))
  s2e <- stats::sigma(model)^2
  C <- if (need_C) {
    rc <- pm_ranef_cov(model)
    i <- which(rc$tab$group == genotype & rc$tab$term == "(Intercept)")
    rc$V[length(rc$beta) + i, length(rc$beta) + i, drop = FALSE]
  }
  usable <- if (is.null(user_vars)) names(dat) else intersect(names(dat), user_vars)
  list(s2g = s2g, terms = s2t, s2e = s2e, n_terms = n_terms, n_plots = n_plots, n_env_arg = n_env, n_rep_arg = n_rep,
       unequal = length(unique(as.vector(kp))) > 1, need_C = need_C, user_vars = user_vars,
       incomplete = pm_h2_incomplete(dat, c(names(vc), pm_h2_fixed_groups(model, dat)), genotype, usable),
       nested_in = pm_h2_nested_fixed(model, dat, gcol), C = C, with_g = with_g, genotype = genotype)
}

## The genotype factor of the model frame; a genotype written as an
## interaction (Family:Province, from (1 | Province/Family)) is made from its parts.
pm_h2_gcol <- function(dat, genotype) {
  if (genotype %in% names(dat)) return(factor(dat[[genotype]]))
  v <- strsplit(genotype, ":", fixed = TRUE)[[1]]
  if (!all(v %in% names(dat))) stop("'", genotype, "' is not a factor of the model.", call. = FALSE)
  interaction(dat[v], drop = TRUE, sep = ":")
}

## Fixed factor terms of the model (blocks fitted as fixed effects).
pm_h2_fixed_groups <- function(model, dat) {
  tl <- attr(stats::terms(pm_nobars(stats::formula(model))), "term.labels")
  tl[vapply(strsplit(tl, ":", fixed = TRUE), function(v) all(v %in% names(dat)) &&
              all(vapply(dat[v], function(x) is.factor(x) || is.character(x), logical(1))), logical(1))]
}

## A fixed factor that the genotypes are nested in (provenances, each family
## in one provenance): its name, or NULL.
pm_h2_nested_fixed <- function(model, dat, gcol) {
  fx <- pm_h2_fixed_groups(model, dat)
  fx <- fx[!grepl(":", fx, fixed = TRUE)]
  for (f in fx) {
    lv <- dat[[f]]
    if (length(unique(lv)) > 1 && all(tapply(lv, gcol, function(x) length(unique(x))) <= 1, na.rm = TRUE)) return(f)
  }
  NULL
}

## Incomplete blocks: a random term without the genotype whose levels hold
## fewer genotypes than the trial (an alpha or lattice design). The standard
## form then ignores the block adjustment of the genotype means.
pm_h2_incomplete <- function(dat, groups, genotype, usable = names(dat)) {
  gcol <- pm_h2_gcol(dat, genotype)
  ng <- length(unique(gcol))
  gv <- strsplit(genotype, ":", fixed = TRUE)[[1]]
  any(vapply(unique(groups), function(g) {
    v <- strsplit(g, ":", fixed = TRUE)[[1]]
    if (any(gv %in% v) || !all(v %in% usable)) return(FALSE)
    k <- tapply(gcol, do.call(paste, c(dat[v], sep = "\r")), function(x) length(unique(x)))
    any(k < ng)
  }, logical(1)))
}

## Mean variance of a difference between the m entries with covariance C:
## mean over pairs of C_ii + C_jj - 2 C_ij = 2 (m tr(C) - sum(C)) / (m (m - 1)).
pm_h2_vdiff <- function(C) {
  m <- nrow(C)
  2 * (m * sum(diag(C)) - sum(C)) / (m * (m - 1))
}

## vBLUE with the variance components fixed at those of the random-genotype
## fit (Schmidt et al. 2019): genotype as a fixed effect, the other random
## terms with their estimated variances, var(beta) = (X'V^-1X)^-1 by the
## Woodbury form V^-1 = W/s2 - W A (I + A'WA)^-1 A'W/s2, A = Z_o Lambda_o
## (sparse; no refit, no n x n matrix). The genotype coefficients are
## differences from the first genotype (treatment coding), which is added
## with variance 0.
pm_h2_vblue_fixed <- function(model, genotype) {
  mf <- stats::model.frame(model)
  mf$.pm_g <- pm_h2_gcol(mf, genotype)
  fixed_terms <- attr(stats::terms(pm_nobars(stats::formula(model))), "term.labels")
  intercept <- attr(stats::terms(pm_nobars(stats::formula(model))), "intercept") == 1
  ff <- stats::reformulate(c(fixed_terms, ".pm_g"), intercept = intercept)
  X0 <- stats::model.matrix(ff, mf, contrasts.arg = list(.pm_g = "contr.treatment"))
  q <- qr(X0)
  X <- X0[, q$pivot[seq_len(q$rank)], drop = FALSE]
  w <- model@resp$weights
  ztl <- lme4::getME(model, "Ztlist")
  rows <- rep(seq_along(ztl), vapply(ztl, nrow, 1L))
  keep <- which(rows != which(names(ztl) == paste0(genotype, ".(Intercept)")))
  s2 <- stats::sigma(model)^2
  if (length(keep)) {
    A <- Matrix::t(lme4::getME(model, "Lambdat")[keep, keep, drop = FALSE] %*% lme4::getME(model, "Zt")[keep, , drop = FALSE])
    WA <- A * w
    M <- Matrix::forceSymmetric(Matrix::crossprod(A, WA) + Matrix::Diagonal(length(keep)))
    XWA <- as.matrix(Matrix::crossprod(X, WA))
    info <- (crossprod(X, X * w) - XWA %*% as.matrix(Matrix::solve(M, Matrix::t(Matrix::Matrix(XWA))))) / s2
  } else {
    info <- crossprod(X, X * w) / s2
  }
  V <- solve(info)
  levs <- levels(mf$.pm_g)
  cols <- paste0(".pm_g", levs)
  # a genotype column left out by the coding (the first level) has effect 0;
  # one dropped as aliased (genotypes nested in a fixed factor) has no BLUE
  if (any(cols %in% colnames(X0) & !cols %in% colnames(X))) return(NA_real_)
  have <- cols %in% colnames(V)
  C <- matrix(0, length(levs), length(levs))
  C[have, have] <- V[cols[have], cols[have]]
  pm_h2_vdiff(C)
}

## vBLUE as in Piepho and Möhring (2007): the same model with genotype fixed
## refitted, its variance components estimated again (its other random terms
## kept), the genotype means (Kmatrix()) and their unadjusted covariance; lm()
## when no random term is left.
pm_h2_vblue <- function(model, genotype, y = NULL) {
  f <- stats::formula(model)
  bars <- pm_findbars(f)
  keep <- bars[vapply(bars, function(b) !identical(deparse(b), paste("1 |", genotype)), logical(1))]
  fixed_terms <- attr(stats::terms(pm_nobars(f)), "term.labels")
  rhs <- c(fixed_terms, ".pm_g", vapply(keep, function(b) paste0("(", deparse(b), ")"), ""))
  # the response as a column (also a transformed one, log(y), or a bootstrap
  # draw); the data of the call (covariates such as log(x) are evaluated
  # there), on the rows the model used; its weights and offset kept
  newf <- stats::as.formula(paste(".pm_y ~", paste(rhs, collapse = " + ")), env = environment(f))
  mf <- stats::model.frame(model)
  dat <- tryCatch(pm_get_data(model), error = function(e) NULL)
  dat <- if (!is.null(dat) && all(rownames(mf) %in% rownames(dat))) dat[rownames(mf), , drop = FALSE] else mf
  dat$.pm_y <- if (is.null(y)) stats::model.response(mf) else y
  dat$.pm_g <- pm_h2_gcol(mf, genotype)
  # (do.call() passes the values: lmer() and lm() evaluate weights and offset
  # among the data's columns)
  args <- list(formula = newf, data = dat, weights = stats::model.weights(mf), offset = stats::model.offset(mf))
  fit <- if (length(keep)) {
    suppressMessages(do.call(lme4::lmer, c(args, REML = lme4::getME(model, "is_REML"))))
  } else {
    do.call(stats::lm, args)
  }
  if (any(is.na(pm_coef(fit)[grep("^\\.pm_g", names(pm_coef(fit)))]))) return(NA_real_)
  km <- suppressMessages(suppressWarnings(Kmatrix(fit, ".pm_g")))
  mp <- pm_model_parm(fit)
  K <- pm_K_cols(km$K, mp)
  pm_h2_vdiff(K %*% mp$vcov %*% t(K))
}

## Simulated heritability of Piepho and Möhring (2007, Eqs. 15-17; Schmidt et
## al. 2019, Eq. 41): the mean over draws of the squared correlation between
## the true genotypic values g and their BLUPs, from their joint covariance
## Omega = (s2g I, M; M, M), M = var(g-hat) = s2g I - C (independent genotypes).
## Drawn in batches of columns, from Omega's eigen decomposition (it is
## singular: the genotype mean is confounded with the intercept).
pm_h2_sim <- function(C, s2g, nsim = 10000) {
  m <- nrow(C)
  M <- s2g * diag(m) - as.matrix(C)
  Om <- rbind(cbind(s2g * diag(m), M), cbind(M, M))
  e <- eigen((Om + t(Om))/2, symmetric = TRUE)
  G <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), 2*m)
  r2 <- numeric(0)
  for (k in split(seq_len(nsim), ceiling(seq_len(nsim) / 2000))) {
    W <- G %*% matrix(stats::rnorm(2*m*length(k)), 2*m)
    g <- W[seq_len(m), , drop = FALSE]; gh <- W[m + seq_len(m), , drop = FALSE]
    gc <- sweep(g, 2, colMeans(g)); hc <- sweep(gh, 2, colMeans(gh))
    r2 <- c(r2, colSums(gc * hc)^2 / (colSums(gc^2) * colSums(hc^2)))
  }
  structure(mean(r2), mc_se = stats::sd(r2) / sqrt(nsim))
}

## Oakey et al. (2006): the mean of the non-zero eigenvalues of I - G^-1 C
## (G = s2g I); the zero ones are directions confounded with fixed effects
## (the mean of the genotypes, with the intercept).
pm_h2_oakey <- function(C, s2g) {
  ev <- eigen(diag(nrow(C)) - as.matrix(C)/s2g, symmetric = TRUE, only.values = TRUE)$values
  mean(ev[ev > 1e-8])
}

## The variance among plants within plots and the number of plants per plot,
## for plot-mean data: numbers, or names of columns of the model's data (the
## mean of the variances, the harmonic mean of the numbers). NULL without them.
pm_h2_within <- function(model, within_variance, within_n) {
  if (is.null(within_variance) && is.null(within_n)) return(NULL)
  if (is.null(within_variance) || is.null(within_n)) {
    stop("Give both 'within_variance' and 'within_n' (plot means: the variance among plants within plots and ",
         "the number of plants per plot).", call. = FALSE)
  }
  get <- function(x) {
    if (is.character(x)) {
      dat <- pm_get_data(model)
      if (!x %in% names(dat)) stop("'", x, "' is not a column of the model's data.", call. = FALSE)
      dat[[x]]
    } else x
  }
  v <- get(within_variance); n <- get(within_n)
  v <- v[is.finite(v)]; n <- n[is.finite(n) & n > 0]
  list(s2w = mean(v), n = length(n) / sum(1/n))
}

## Names of random intercept terms as pm_vc_mer_fun() gives them (none for none:
## paste0() would return ".(Intercept)" for an empty vector).
pm_h2_tnames <- function(t) if (length(t)) paste0(t, ".(Intercept)") else character(0)

## SE of a form that is a function of the variance components, by the delta
## method on the variance parameters (theta, sigma) with their covariance from
## the observed information (pm_lmer_varpar(), as varcomp()).
pm_h2_delta <- function(model, f) {
  fun <- pm_vc_mer_fun(model, residual = TRUE)
  vp <- pm_lmer_varpar(model)
  J <- numDeriv::jacobian(function(x) f(fun(x)), vp$varpar)
  sqrt(drop(J %*% vp$vcov_varpar %*% t(J)))
}

## Parametric bootstrap of the forms that depend on prediction errors: draws
## from the fitted model, each refitted (pm_mer_refit()) and its forms
## recomputed. Seeds per draw are made here, so the result is the same for any
## ncore (pm_perm_map(), R/perm-engine.R). Returns a matrix nboot x methods.
pm_h2_boot <- function(model, comp, methods, nboot, ncore = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, nboot)
  ys <- lapply(seeds, function(s) { set.seed(s); stats::simulate(model, nsim = 1)[[1]] })
  one <- function(b) {
    fit <- suppressMessages(suppressWarnings(pm_mer_refit(model, ys[[b]])))
    cb <- tryCatch(pm_h2_components(fit, comp$genotype, comp$n_env_arg, comp$n_rep_arg, need_C = comp$need_C,
                                    user_vars = comp$user_vars), error = function(e) NULL)
    if (!is.null(cb)) { cb$piepho <- comp$piepho; cb$nsim <- comp$nsim }
    if (is.null(cb)) return(stats::setNames(rep(NA_real_, length(methods)), methods))
    set.seed(seeds[b])                                  # the simulated form draws too
    pm_h2_point(fit, cb, methods, y = ys[[b]])
  }
  out <- pm_perm_map(seq_len(nboot), one, ncore, packages = "predictmeans")
  out <- lapply(out, function(r) if (is.null(r)) stats::setNames(rep(NA_real_, length(methods)), methods) else r)
  do.call(rbind, out)
}

## Limits on the logit scale of est x scale, so they stay in (0, 1/scale)
## (scale = relationship for the family form, which ranges up to 1/relationship).
pm_h2_limits <- function(est, se, alpha, scale = 1) {
  e <- est * scale
  if (is.na(se) || is.na(e) || e <= 0 || e >= 1) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - alpha/2)
  l <- stats::qlogis(e) + c(-1, 1) * z * se * scale / (e * (1 - e))
  stats::plogis(l) / scale
}

#' @exportS3Method print pm_heritability
print.pm_heritability <- function(x, digits = 4, ...) {
  cmp <- attr(x, "components")
  if (is.null(cmp)) {                        # a subset ([ drops the attributes): a plain table
    print(structure(x, class = "data.frame"), digits = digits, row.names = FALSE)
    return(invisible(x))
  }
  m <- x$method
  kind <- c(if (any(m %in% c("standard", "cullis", "piepho", "oakey", "reliability", "simulated"))) "heritability on a line-mean basis",
            if ("family" %in% m) paste0("narrow-sense heritability of a single plant (relationship ",
                                        format(attr(x, "relationship"), digits = 3), ")"),
            if ("repeatability" %in% m) "repeatability of a single record")
  cat(toupper(substring(kind[1], 1, 1)), substring(paste(kind, collapse = "; "), 2), " (",
      format(100 * attr(x, "conf_level")), "% limits)\n", sep = "")
  cat("Variance of ", attr(x, "genotype"), " ", format(cmp$s2g, digits = digits), "; residual ",
      format(cmp$s2e, digits = digits),
      if (!is.null(cmp$within) && any(m %in% c("family", "repeatability")))
        paste0("; within plots ", format(cmp$within$s2w, digits = digits), " (", format(cmp$within$n, digits = 3),
               " plants per plot)"),
      if (length(cmp$terms)) paste0("; ", paste(names(cmp$terms), format(cmp$terms, digits = digits), sep = " ",
                                                 collapse = ", ")), "\n\n", sep = "")
  y <- x
  class(y) <- "data.frame"
  attr(y, "components") <- attr(y, "conf_level") <- NULL
  print(y, digits = digits, row.names = FALSE)
  if (!is.null(attr(x, "mc_se"))) {
    cat("\nSimulated form: Monte Carlo SE ", format(attr(x, "mc_se"), digits = 2), " (", attr(x, "nsim"),
        " draws; nsim for more).\n", sep = "")
  }
  nb <- attr(x, "nboot")
  if (isTRUE(nb > 0)) cat("\nSE and limits of the mixed-model forms from", nb, "parametric bootstrap draws (percentile).\n")
  pm_print_advice(x)
  invisible(x)
}
