## Help pages of the data sets (roxygen2 sources).

#' \code{ATP} containing data
#'
#' @description
#' \code{ATP} containing data from an experiment to study the effects 
#' of preserving liquids on the enzyme content of dog hearts. There were 23 hearts 
#' and two treatment factors, A and B, each at two levels. Measurements were made 
#' of ATP as a percentage of total enzyme in the heart, at one and two hourly
#' intervals during a twelve hour period following initial preservation.
#' @usage
#' data(ATP)
#' @format
#' \code{ATP} is a 230 row data frame with the following columns 
#' \describe{
#' \item{heart}{dog heart id.}
#' \item{time}{time in hour for \code{ATP} measurement.}
#' \item{A}{treatment with two levels.}
#' \item{B}{treatment with two levels.}
#' \item{ATP}{percentage of total enzyme in the heart.}
#' }
#' @name ATP
#' @docType data
NULL

#' \code{Body_length} data
#'
#' @description
#' Simulated body lengths of beetles reared in containers within populations, under two
#' treatments in two habitats, of the kind used by Nakagawa and Schielzeth (2013) to illustrate
#' \eqn{R^2} for mixed models. Used in the \code{R2_glmm} vignette.
#' @usage
#' data(Body_length)
#' @format
#' A data frame with 960 rows and 6 columns:
#' \describe{
#' \item{Population}{factor, 12 populations.}
#' \item{Container}{factor, 120 containers (nested in populations).}
#' \item{Sex}{factor, \code{Female} or \code{Male}.}
#' \item{Habitat}{factor, \code{A} or \code{B}.}
#' \item{Treatment}{factor, \code{Cont} or \code{Exp}.}
#' \item{BodyL}{numeric, body length.}
#' }
#' @references
#' Nakagawa, S. and Schielzeth, H. (2013), A general and simple method for obtaining
#' \eqn{R^2} from generalized linear mixed-effects models. \emph{Methods in Ecology and
#' Evolution}, 4, 133--142.
#' @examples
#' data(Body_length)
#' fm <- lmer(BodyL ~ Treatment + Habitat + Sex + (1|Population/Container), data = Body_length)
#' R2_glmm(fm)
#' @name Body_length
#' @docType data
NULL

#' \code{Clinical} data
#'
#' @description
#' \code{Clinical} data is from a multicentre randomized clinical trial (Beitler
#'  & Landis 1985, Biometrics).
#' @usage
#' data(Clinical)
#' @format
#' \code{Clinical} is a 16 row data frame with the following columns 
#' \describe{
#' \item{Clinic}{8 centres id.}
#' \item{Treatment}{2 skin treatments (control or active drug).}
#' \item{Favorable}{number that produced a favourable response.}
#' \item{Total}{number of patients in each treatment group.}
#' }
#' @name Clinical
#' @docType data
NULL

#' \code{Drug} data
#'
#' @description
#' The data is for the comparison of the effectiveness of three analgesic drugs 
#' to a standard drug, morphine (Finney, Probit analysis, 3rd Edition 1971, p.103). 14 groups 
#' of mice were tested for response to the drugs at a range of doses.
#' @usage
#' data(Drug)
#' @format
#' \code{Drug} is a 14 row data frame with the following columns 
#' \describe{
#' \item{Drug}{type of drug.}
#' \item{Dose}{dose volumn.}
#' \item{N}{total number of mice in each group.}
#' \item{R}{number responding mice in each group.}
#' \item{log10Dose}{log10 transformed dose volumn.}
#' }
#' @name Drug
#' @docType data
NULL

#' \code{Poisson_ind} data
#'
#' @description
#' Simulated egg counts of individual female beetles in containers within populations, under
#' two treatments in two habitats, of the kind used by Nakagawa and Schielzeth (2013). Used
#' in the \code{R2_glmm} vignette for a Poisson GLMM with an individual-level (overdispersion)
#' random effect.
#' @usage
#' data(Poisson_ind)
#' @format
#' A data frame with 480 rows and 7 columns:
#' \describe{
#' \item{Population}{factor, 12 populations.}
#' \item{Container}{factor, 60 containers (nested in populations).}
#' \item{Sex}{factor, \code{Female} only.}
#' \item{Habitat}{factor, \code{B} or \code{A}.}
#' \item{Treatment}{factor, \code{Exp} or \code{Cont}.}
#' \item{Egg}{integer, number of eggs.}
#' \item{individual}{factor, individual within container.}
#' }
#' @references
#' Nakagawa, S. and Schielzeth, H. (2013), A general and simple method for obtaining
#' \eqn{R^2} from generalized linear mixed-effects models. \emph{Methods in Ecology and
#' Evolution}, 4, 133--142.
#' @name Poisson_ind
#' @docType data
NULL

#' \code{copper} data
#'
#' @description
#' Copper grade at sampling points with their X, Y and Z coordinates within a designated
#' area. Used in the overview article to fit a semiparametric model with a smooth in Z and a
#' thin-plate spline in X and Y.
#' @usage
#' data(copper)
#' @format
#' A data frame with 439 rows and 4 columns:
#' \describe{
#' \item{xcoord}{X coordinate of the sampling point.}
#' \item{ycoord}{Y coordinate of the sampling point.}
#' \item{zcoord}{Z coordinate of the sampling point.}
#' \item{grade}{copper grade.}
#' }
#' @examples
#' data(copper)
#' \donttest{
#' mod <- semireg(log(grade) ~ (xcoord + ycoord + zcoord)^2,
#'                smoothZ = list(grp_z = smZ(zcoord),
#'                               grp_geo = smZ(cbind(xcoord, ycoord), type = "Ztps")),
#'                data = copper, prt = FALSE)
#' ranova(mod$semer)
#' }
#' @name copper
#' @docType data
NULL
