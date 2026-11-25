#' Standard Errors of the Predicted Random Effects
#'
#' @description
#' The predicted random effects (BLUPs, conditional modes) of a mixed model with their standard errors and
#' limits, and optionally the predicted value of each level (fixed part plus random effect).
#'
#' Two standard errors are available:
#' \describe{
#' \item{\code{type = "prediction"}}{(the default) the prediction error SE, \eqn{\sqrt{Var(\hat u - u)}}, from
#'   the inverse of the mixed model equations (Henderson's \eqn{C_{22}}; SAS PROC MIXED \code{RANDOM / SOLUTION},
#'   "Std Err Pred"). It includes the uncertainty of the estimated fixed effects, which matters most for a term
#'   that contains a fixed factor (e.g. \code{Block:Variety} with \code{Variety} fixed).}
#' \item{\code{type = "conditional"}}{the conditional SE, \eqn{\sqrt{Var(u \mid y)}} with the fixed effects
#'   treated as known: lme4's \code{ranef(condVar = TRUE)}, the numbers of \code{se_ranef()} in predictmeans 1.x.
#'   It is smaller.}
#' }
#' Both treat the variance parameters as known (the usual plug-in), except for \code{glmmTMB}: there
#' \code{"conditional"} is \code{ranef(condVar = TRUE)} of glmmTMB, which is TMB's \code{sdreport()} and
#' includes the uncertainty of the fixed effects \emph{and} of the variance parameters (so it can exceed
#' \code{"prediction"}, which inverts the \eqn{(\beta, b)} block of TMB's joint precision, variance parameters
#' known). Correlated terms, nested and crossed terms, glmmTMB structures (\code{us}, \code{cs}, \code{ar1},
#' ...) and \code{lme} correlation and variance structures are all covered, as the SEs come from the joint
#' covariance of all random effects. For GLMMs, all values are on the link scale.
#'
#' With \code{predict = TRUE}, each level gets its predicted value: for a random term with no fixed factor in it
#' (\code{Block}), the overall mean (the fixed factors averaged with equal weights, as \code{\link{predictmeans}};
#' covariates at their means) plus \eqn{\hat u}; for a term that contains fixed factors (\code{Block:Variety},
#' \code{Variety} fixed), the predicted mean of that fixed level plus \eqn{\hat u}. Its SE is from the joint
#' covariance of \eqn{(\hat\beta, \hat u - u)}, \eqn{L C L'} with \eqn{L = [k', e_i]} (\code{"prediction"}
#' always). Only intercept rows get a prediction (random slopes: \code{NA}). For GLMMs the response scale is added
#' by the inverse link (SE by the delta method).
#' @param object An \code{lmer}, \code{glmer}, \code{glmmTMB} (conditional model), \code{lme} or
#'   \code{\link{semireg}} fit.
#' @param rand_term Optional name of one grouping factor, e.g. \code{"Block"} or \code{"Variety:Block"}
#'   (\code{lme}: the level name, e.g. \code{"Variety"} for \code{~ 1 | Block/Variety}).
#' @param type \code{"prediction"} (default) or \code{"conditional"} (1.x); see Description.
#' @param predict Add the predicted value of each level (see Description)?
#' @param conf_level Confidence level of the limits.
#' @param x A \code{pm_ranef} result.
#' @param n Rows printed per group.
#' @param digits Digits printed.
#' @param ... Not used.
#' @return
#' A data frame of class \code{pm_ranef}: \code{group}, \code{level}, \code{term}, \code{estimate} (the BLUP),
#' \code{std.error}, \code{conf.low}, \code{conf.high} (normal limits \eqn{\hat u \pm z\,SE}: for
#' \code{"prediction"} a prediction interval of \eqn{u}); with \code{predict = TRUE} also \code{prediction},
#' \code{prediction.std.error}, \code{prediction.conf.low}, \code{prediction.conf.high}, and for GLMMs
#' \code{response}, \code{response.std.error}, \code{response.conf.low}, \code{response.conf.high}.
#' \code{attr(, "settings")} records \code{type}, \code{conf_level}, \code{limits}, \code{model_class},
#' \code{scale} and \code{conditional} (what the conditional SE is for this class). \code{as.list()} gives the 1.x
#' result: one matrix of SEs (levels x terms) per group, or that matrix only when \code{rand_term} was given.
#' @references
#' Henderson, C. R. (1975). Best linear unbiased estimation and prediction under a selection model.
#' \emph{Biometrics}, 31, 423-447.
#' @seealso \code{\link{varcomp}}, \code{\link{residplot}}, \code{\link{advice}}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#' se_ranef(fm, "Block")
#' se_ranef(fm, "Block", type = "conditional")     # 1.x numbers
#' se_ranef(fm, "Variety:Block", predict = TRUE)   # cell values
#' as.list(se_ranef(fm, type = "conditional"))     # 1.x list of matrices
#' @aliases print.pm_ranef tidy.pm_ranef as.data.frame.pm_ranef as.list.pm_ranef
#' @export
se_ranef <- function(object, rand_term = NULL, type = c("prediction", "conditional"), predict = FALSE,
                     conf_level = 0.95) {
  type <- match.arg(type)
  stopifnot(is.numeric(conf_level), length(conf_level) == 1, conf_level > 0, conf_level < 1)
  rc <- pm_ranef_cov(object)
  tab <- rc$tab
  groups <- unique(tab$group)
  if (!is.null(rand_term)) {
    if (!is.character(rand_term) || length(rand_term) != 1 || !rand_term %in% groups) {
      stop("'rand_term' must be one of the grouping factors: ", paste0("\"", groups, "\"", collapse = ", "), ".",
           call. = FALSE)
    }
  }
  p <- length(rc$beta)
  V <- rc$V
  q <- length(rc$u)
  vu <- if (type == "prediction") diag(V)[p + seq_len(q)] else rc$cond
  z <- stats::qnorm(1 - (1 - conf_level)/2)
  out <- data.frame(group = tab$group, level = tab$level, term = tab$term, estimate = rc$u,
                    std.error = sqrt(pmax(vu, 0)), stringsAsFactors = FALSE)
  out$conf.low <- out$estimate - z*out$std.error
  out$conf.high <- out$estimate + z*out$std.error
  keep <- if (is.null(rand_term)) rep(TRUE, q) else tab$group == rand_term
  link <- rc$family
  resp <- !is.null(link) && !identical(link$link, "identity")
  if (predict) out <- pm_ranef_predict(object, rc, out, keep, z, resp)
  out <- out[keep, , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("pm_ranef", "data.frame")
  attr(out, "settings") <- list(
    type = type, conf_level = conf_level, rand_term = rand_term, predict = predict,
    limits = paste0("normal ", if (type == "prediction") "prediction intervals of u" else
      "intervals from the conditional SE"),
    model_class = class(pm_fit(object))[1],
    scale = if (!resp) "identity" else paste0("link (", link$link, ")"),
    conditional = rc$cond_label)
  pm_add_advice(out, pm_ranef_advice(object, rc, out, type))
}

## The predicted value of each level (intercept rows): k'beta + u, SE from the
## joint covariance with L = [k', e_i].
pm_ranef_predict <- function(object, rc, out, keep, z, resp) {
  p <- length(rc$beta)
  V <- rc$V
  out$prediction <- out$prediction.std.error <- NA_real_
  for (g in unique(out$group[keep])) {
    rows_g <- which(out$group == g & out$term == "(Intercept)")
    if (!length(rows_g)) next
    gf <- rc$rows[[g]]
    vars <- intersect(all.vars(str2lang(paste0("~", gsub("/", ":", g)))), names(rc$data))
    if (inherits(pm_fit(object), "lme")) {
      vars <- intersect(all.vars(nlme::getGroupsFormula(object)), names(rc$data))
      vars <- vars[seq_len(match(g, vars))]
    }
    fx <- intersect(vars, rc$fixed)
    if (length(fx)) {
      km <- Kmatrix(object, paste(fx, collapse = ":"))
      K <- km$K
      key <- do.call(paste, c(lapply(fx, function(v) as.character(rc$data[[v]])), sep = ":"))
      key <- key[match(out$level[rows_g], as.character(gf))]
      rn <- do.call(paste, c(lapply(fx, function(v) as.character(km$fctnames[[v]])), sep = ":"))
      Kr <- K[match(key, rn), names(rc$beta), drop = FALSE]
    } else {
      K <- if (length(rc$fixed)) colMeans(Kmatrix(object, rc$fixed[1])$K[, names(rc$beta), drop = FALSE]) else {
        X <- if (inherits(object, "lme")) stats::model.matrix(stats::formula(object), rc$data) else
          lme4::getME(pm_fit(object), "X")
        colMeans(as.matrix(X))[names(rc$beta)]
      }
      Kr <- matrix(K, length(rows_g), p, byrow = TRUE)
    }
    ui <- p + rows_g
    V11 <- V[seq_len(p), seq_len(p)]
    v <- rowSums((Kr %*% V11) * Kr) + 2*rowSums(Kr * t(V[seq_len(p), ui, drop = FALSE])) + diag(V)[ui]
    out$prediction[rows_g] <- as.vector(Kr %*% rc$beta) + out$estimate[rows_g]
    out$prediction.std.error[rows_g] <- sqrt(pmax(v, 0))
  }
  out <- out[, c(setdiff(names(out), c("prediction", "prediction.std.error")), "prediction", "prediction.std.error")]
  out$prediction.conf.low <- out$prediction - z*out$prediction.std.error
  out$prediction.conf.high <- out$prediction + z*out$prediction.std.error
  if (resp) {
    f <- rc$family
    out$response <- f$linkinv(out$prediction)
    out$response.std.error <- abs(f$mu.eta(out$prediction))*out$prediction.std.error
    lo <- f$linkinv(out$prediction.conf.low)
    hi <- f$linkinv(out$prediction.conf.high)
    out$response.conf.low <- pmin(lo, hi)
    out$response.conf.high <- pmax(lo, hi)
  }
  out
}

## Advice: a variance at zero, few levels, the conditional SE of a term that
## contains a fixed factor.
pm_ranef_advice <- function(object, rc, out, type) {
  rows <- list()
  tab <- rc$tab[rc$tab$group %in% out$group, , drop = FALSE]
  vz <- unique(tab$group[tab$variance < 1e-8 * max(1, abs(rc$tab$variance))])
  if (length(vz)) {
    rows$zero <- pm_advice("variance_at_zero", paste0(
      "A variance of ", paste(vz, collapse = ", "), " is estimated as zero: its predicted effects are all 0 with SE 0, ",
      "which does not show that the levels are equal. Check it with varcomp(test = TRUE)."), "caution")
  }
  nl <- vapply(split(tab$level, tab$group), function(l) length(unique(l)), 1L)
  few <- nl[nl < 5]
  if (length(few)) {
    rows$few <- pm_advice("few_re_levels", paste0(
      "Few levels for ", paste0(names(few), " (", few, ")", collapse = ", "), ": its variance, and so the shrinkage ",
      "of its predicted effects, is poorly estimated."), "note")
  }
  if (type == "conditional") {
    conf <- unique(out$group)[vapply(unique(out$group), function(g) {
      length(intersect(all.vars(str2lang(paste0("~", gsub("/", ":", g)))), rc$fixed)) > 0 ||
        (inherits(object, "lme") && length(intersect(all.vars(nlme::getGroupsFormula(object))[
          seq_len(match(g, all.vars(nlme::getGroupsFormula(object)), nomatch = 0))], rc$fixed)) > 0)
    }, logical(1))]
    if (length(conf)) {
      rows$cond <- pm_advice("ranef_conditional_fixed", paste0(
        paste(conf, collapse = ", "), " contains a fixed factor: the conditional SE treats its fixed effects as ",
        "known and understates the error of the predicted effects; type = \"prediction\" includes it."), "caution")
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}

#' @rdname se_ranef
#' @exportS3Method print pm_ranef
print.pm_ranef <- function(x, n = 10, digits = 4, ...) {
  s <- attr(x, "settings")
  cat("Predicted random effects (BLUPs), ", if (s$type == "prediction") "prediction error SE" else
    "conditional SE", "; ", format(100*s$conf_level), "% ", s$limits, if (s$scale != "identity") paste0("; ", s$scale, " scale"), "\n", sep = "")
  tab <- as.data.frame(x)
  for (g in unique(tab$group)) {
    tg <- tab[tab$group == g, setdiff(names(tab), "group"), drop = FALSE]
    cat("\n", g, " (", length(unique(tg$level)), " levels)\n", sep = "")
    print(utils::head(tg, n), digits = digits, row.names = FALSE)
    if (nrow(tg) > n) cat("... ", nrow(tg) - n, " more rows: as.data.frame(x) gives the table.\n", sep = "")
  }
  pm_print_advice(x)
  invisible(x)
}

#' @exportS3Method as.data.frame pm_ranef
as.data.frame.pm_ranef <- function(x, row.names = NULL, optional = FALSE, ...) {
  attr(x, "settings") <- attr(x, "advice") <- NULL
  class(x) <- "data.frame"
  if (!is.null(row.names)) rownames(x) <- row.names
  x
}

#' @exportS3Method generics::tidy pm_ranef
tidy.pm_ranef <- function(x, ...) as.data.frame.pm_ranef(x)

#' @exportS3Method as.list pm_ranef
as.list.pm_ranef <- function(x, ...) {
  tab <- as.data.frame.pm_ranef(x)
  out <- lapply(split(tab, factor(tab$group, unique(tab$group))), function(d) {
    lev <- unique(d$level)
    tm <- unique(d$term)
    m <- matrix(NA_real_, length(lev), length(tm), dimnames = list(lev, tm))
    m[cbind(match(d$level, lev), match(d$term, tm))] <- d$std.error
    m
  })
  if (!is.null(attr(x, "settings")$rand_term)) out[[1]] else out
}
