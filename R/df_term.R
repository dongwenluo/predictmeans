#' Calculate degree of freedom of a modelterm (contrast) for a lmer model
#'
#' Calculate the degree of freedom of a modelterm (contrast) for a \code{lmer}
#' model using "Kenward-Roger" or "Satterthwaite" method.
#'
#'
#' @param model Model object returned by \code{lmer}.
#' @param modelterm Name (in "quotes") for indicating which factor term's
#' degree of freedom to be calculated.  The \code{modelterm} must be given
#' exactly as it appears in the model formlar, e.g. "A" or "A:B".
#' @param covariate Name (in "quotes") of one the covariate variables in the
#' \code{model}.
#' @param ctrmatrix A specified contrast matrix. If \code{ctrmatrix} isn't
#' NULL, the programe will ignore modelterm and calculate degree of freedom for
#' the \code{ctrmatrix}.
#' @param ctrnames Names of the specified contrasts, e.g. c("A vs D", "C vs B",
#' ...)
#' @param type Name (in "quote") for indicating a method for claculating degree
#' of freedom.  The choices are "Kenward-Roger" and "Satterthwaite". The
#' default method is "Kenward-Roger".
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'
#' library(predictmeans)
#' # ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#' df_term(fm, "nitro:Variety")
#' #--------------------------------------------------------
#' ## Not run:
#' ## The contrast has a contrast matrix as follows:
#' #     0:Golden Rain 0:Marvellous 0:Victory
#' #[1,]            -1            0         1
#' #[2,]             0            0         1
#' #     0.2:Golden Rain 0.2:Marvellous 0.2:Victory
#' #[1,]               0              0           0
#' #[2,]               0              0           0
#' #     0.4:Golden Rain  0.4:Marvellous 0.4:Victory
#' #[1,]               0               0           0
#' #[2,]               0              -1           0
#' #      0.6:Golden Rain 0.6:Marvellous 0.6:Victory
#' #[1,]                0              0           0
#' #[2,]                0              0           0
#' #--------------------------------------------------------
#' # 1. Enter above contrast matrix into a pop up window, then close the window
#' # df_term(fm, "nitro:Variety")
#'
#' # 2. Construct the contrast matrix directly
#' cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'             c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#' df_term(fm, ctrmatrix=cm, type="Satterthwaite")
#'
#' @importFrom lmerTest as_lmerModLmerTest
#' @export

df_term <- function(model,
                    modelterm,
                    covariate = NULL,
                    ctrmatrix = NULL,
                    ctrnames = NULL,
                    type = c("Kenward-Roger", "Satterthwaite")) {
  stopifnot(inherits(model, "lmerMod"))

  if (!getME(model, "is_REML")) {
    stop("This function works properly only for REML model fits")
  }
  if (!is.null(ctrmatrix)) {
    if (is.vector(ctrmatrix)) {
      Lc <- matrix(ctrmatrix, nrow = 1)
    } else {
      Lc <- ctrmatrix
    }
    stopifnot(is.numeric(Lc), ncol(Lc) == length(fixef(model)))

    if (!is.null(ctrnames)) {
      rownames(Lc) <- ctrnames
    }
  } else {
    Lc <- Kmatrix(model, modelterm, covariate)$K
  }

  type <- as.character(type)
  type <- match.arg(type)

  if (type == "Kenward-Roger") {
    vcov_beta_adj <- try(pbkrtest::vcovAdj(model), silent = TRUE)
    # Adjusted vcov(beta)
    ddf <- try(apply(Lc, 1, function(x) {
      pbkrtest::Lb_ddf(x, V0 = vcov(model), Vadj = vcov_beta_adj)
    }), silent = TRUE)
    # vcov_beta_adj need to be dgeMatrix!

    if (any(
      inherits(vcov_beta_adj, "try-error"),
      inherits(ddf, "try-error"),
      ddf >= nrow(model.frame(model)),
      ddf <= 0
    )) {
      warning("Unable to compute Kenward-Roger Df: using Satterthwaite instead")
      type <- "Satterthwaite"
    }
  }

  if (type == "Satterthwaite") {
    if (!inherits(model, "lmerModLmerTest")) {
      model <- as_lmerModLmerTest(model)
    }
    ddf <- apply(Lc, 1, function(x) {
      suppressMessages(lmerTest::calcSatterth(model, x)$denom)
    })
  }
  return(ddf)
}
