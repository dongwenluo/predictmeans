# Help page df_term (roxygen2 source; man/df_term.Rd is generated).
#' Calculate degree of freedom of a modelterm (contrast) for a lmer model
#'
#' @description
#' Calculate the degree of freedom of a modelterm (contrast) for a \code{lmer} model 
#'   using "Kenward-Roger" or "Satterthwaite" method.
#' @usage
#' df_term(model, modelterm, covariate=NULL, ctrmatrix=NULL, ctrnames=NULL, 
#'           type=c("Kenward-Roger", "Satterthwaite", "containment"), vcov_adj=NULL)
#' @param model Model object returned by \code{lmer}.
#' @param modelterm Name (in "quotes") for indicating which factor term's degree of freedom to be calculated. 
#'   The \code{modelterm} must be given exactly as it appears in the model formlar, e.g. "A" or "A:B".
#' @param covariate Name (in "quotes") of one the covariate variables in the \code{model}.
#' @param ctrmatrix A specified contrast matrix. If \code{ctrmatrix} isn't NULL, the programe will ignore modelterm 
#'   and calculate degree of freedom for the \code{ctrmatrix}.
#' @param ctrnames Names of the specified contrasts, e.g. c("A vs D", "C vs B", ...)
#' @param type Name (in "quote") for indicating a method for claculating degree of freedom. 
#'   The choices are "Kenward-Roger", "Satterthwaite" and "containment". The default method is "Kenward-Roger".
#'   Kenward-Roger needs a REML fit; for models fitted by ML, Satterthwaite is used with a warning.
#'   "containment" gives the SAS \code{DDFM=CONTAIN} df of \code{modelterm} (needed then) for every row.
#' @param vcov_adj Optional Kenward-Roger adjusted covariance matrix of the fixed effects, as returned by
#'   \code{pbkrtest::vcovAdj(model)}, to avoid computing it again when \code{df_term} is called several times
#'   for the same model.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' # ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#' df_term(fm, "nitro:Variety")
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
#'
#' # 1. Enter above contrast matrix into a pop up window, then close the window
#' if (interactive()) df_term(fm, "nitro:Variety")
#'  
#' # 2. Construct the contrast matrix directly
#' cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#'             c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#' df_term(fm, ctrmatrix=cm, type="Satterthwaite")
#' @name df_term
NULL
