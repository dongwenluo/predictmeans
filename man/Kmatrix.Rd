\name{Kmatrix}
\alias{Kmatrix}
\title{
Matrix of Coefficients in a Linear Model
}
\description{
  This function obtains a matrix of coefficients (the linear functions \eqn{K} such that the
  predicted means are \eqn{K\hat\beta}) for models such as \code{aov}, \code{lm},
  \code{glm}, \code{gls}, \code{lme}, \code{lmer}, \code{glmer}, \code{glmmTMB}, \code{gam}
  (package \pkg{mgcv}) and \code{semireg}.

  The predicted means are built on a reference grid: every combination of the factor levels in
  the model, with covariates fixed at their means (or other chosen values). Each row of \code{K}
  is the equally weighted average of the grid rows for one level combination of
  \code{modelterm}. The design matrix of the grid is evaluated with the model's own stored basis
  information, so covariates may enter the model through transformations such as
  \code{poly(x, 2)}, \code{ns(x, 3)}, \code{bs(x)}, \code{log(x)}, \code{scale(x)} or
  \code{I(x^2)}, only through interactions (e.g. \code{A:x} without \code{x}, or \code{x:z}),
  or through smooth terms (\code{s(x)}, \code{te(x, z)}, \code{s(x, by = A)}) of a \code{gam}.
}
\usage{
Kmatrix(model, modelterm, covariate=NULL, covariateV=NULL, data=NULL, prtnum=FALSE,
        cov_reduce=mean, exclude=NULL)
}
\arguments{
  \item{model}{
Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, \code{lmer},
\code{glmer}, \code{glmmTMB}, \code{gam}/\code{bam} (package \pkg{mgcv}; for \code{gamm} use its
\code{$gam} component) or \code{semireg}.
}
  \item{modelterm}{
  Name (in "quotes") for indicating which model term's predicted mean to be calculated.
  The \code{modelterm} must be given exactly as it appears in the printed model, e.g. "A" or "A:B".
  A factor created in the model formula, e.g. \code{factor(dose)}, can be given either as
  "factor(dose)" or as "dose".
}
  \item{covariate}{
A numerical vector to specify values of covariates for calculating predicted means,
default values are the means of the associated covariates. It also can be the name of one or more covariates in the model,
in which case predictions are made over a range of values of those covariates (see \code{covariateV}).
For a \code{semireg} model this may also be a covariate that only appears in the smooth terms.
}
  \item{covariateV}{
  Values for the covariate(s) named in \code{covariate}. For one covariate, a numeric vector. For several covariates,
  either a list of numeric vectors (all combinations are used) or a matrix with one named column per covariate,
  where each row is one prediction point (e.g. a set of locations for a two-dimensional smooth).
  When \code{covariateV} is \code{NULL}, 50 (one covariate) or 30 (several covariates) equally spaced values over the observed range are used.
}
  \item{data}{
In some cases, you need to provide the data set used in model fitting, for example when the data
used to fit the model can no longer be found. Usually the data are recovered from the model.}
  \item{prtnum}{
An option for printing covariate info on the screen or not. The default is FALSE.
}
  \item{cov_reduce}{
How covariates that are not varied are handled. A function returning a single value, applied to each covariate
(the default \code{mean}; \code{median} is another choice), or the string \code{"average"}, which averages the
predictions over the observed values of those covariates instead of fixing them at one value. With
\code{"average"}, nonlinear terms such as \code{poly(x, 2)} or \code{s(x)} are averaged over the data, which for a
\code{gam} gives predicted means adjusted for (rather than evaluated at the centre of) a spatial or temporal trend.
}
  \item{exclude}{
For \code{gam} models only: labels of smooth terms to leave out of the predictions (their columns of \code{K} are
set to zero), e.g. \code{"s(Block)"}. The default \code{NULL} excludes random-effect smooths (\code{bs = "re"} and
\code{bs = "fs"}), so the predicted means are population-level means, as for \code{lmer}. Use
\code{exclude = character(0)} to keep all smooth terms.
}
}

\value{
 \item{K}{Coefficients matrix. It has an attribute \code{"estimable"}, a logical vector with one value per row that is
 \code{FALSE} for rows that are not estimable functions of the model parameters (e.g. means for a missing cell in a
 rank-deficient model).}
 \item{fctnames}{A model frame contains factor(s) info in the model.}
 \item{response}{The name of response variable in the model.}
 \item{preddf}{The data of the predictor variables used in model fitting.}
 \item{refgrid}{The reference grid used to build \code{K}: a list with the grid data (\code{grid}), the sparse
 weight matrix \code{W} that averages grid rows into rows of \code{K}, and the rows of \code{fctnames}
 (\code{target}). Any other basis evaluated on \code{grid}, e.g. the smooth terms of a \code{semireg} model, can be averaged
 in the same way as \code{W \%*\% basis}.}
}
\references{
  Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004),
	\emph{Prediction in linear mixed models}, Australian and New Zealand Journal of Statistics, 46(3), 325-347.
}
\author{
This function heavily depends on the codes from package "lsmeans".
}

\examples{
  library(predictmeans)
  data(Oats, package="nlme")
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
  Kmatrix(fm, "Variety", prtnum=TRUE)$K
  Kmatrix(fm, "Variety", 0.5, prtnum=TRUE)$K
 # Kmatrix(fm, "Variety", "nitro")$K
  Kmatrix(fm, "Variety", "nitro", covariateV=seq(0, 0.6, 0.1))$K

  # Covariates entering through transformations or only through interactions
  fm2 <- lm(yield ~ Block + Variety + poly(nitro, 2), data=Oats)
  Kmatrix(fm2, "Variety")$K
  fm3 <- lm(yield ~ Variety + Variety:nitro, data=Oats)
  Kmatrix(fm3, "Variety", "nitro", covariateV=c(0, 0.3, 0.6))$K

  # A smooth term in an mgcv model, averaged over the observed nitro values
  if (requireNamespace("mgcv", quietly = TRUE)) {
    fm4 <- mgcv::gam(yield ~ Variety + s(nitro, k = 4) + s(Block, bs = "re"), data = Oats)
    Kmatrix(fm4, "Variety", cov_reduce = "average")$K
  }
}
