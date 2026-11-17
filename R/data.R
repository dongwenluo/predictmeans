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
#' data(copperKnots)
#' @format
#' A data frame with 439 rows and 4 columns:
#' \describe{
#' \item{xcoord}{X coordinate of the sampling point.}
#' \item{ycoord}{Y coordinate of the sampling point.}
#' \item{zcoord}{Z coordinate of the sampling point.}
#' \item{grade}{copper grade.}
#' }
#' \code{copperKnots}: the 20 knots (columns \code{xcoord}, \code{ycoord}) of M. P. Wand's thin-plate spline
#' fit of these data (a geoadditive model fitted by \code{lme()}), which \code{\link{semireg}} reproduces.
#' @examples
#' data(copper)
#' data(copperKnots)
#' # M. P. Wand's geoadditive model: a thin-plate surface over the plane with the given knots,
#' # and an O'Sullivan spline in depth (their linear terms are added to the fixed effects)
#' mod <- semireg(grade ~ s(xcoord, ycoord, bs = "ztps", xt = list(knots = as.matrix(copperKnots))) +
#'                  s(zcoord, bs = "os", xt = list(knots = seq(80, 120, by = 5))), data = copper, prt = FALSE)
#' anova(mod)   # fixed effects (Kenward-Roger) and a test of each smooth term
#' semipred(mod, covariate = c("xcoord", "ycoord"), prt = FALSE)$plt
#' \donttest{
#' # on the log scale, with the pairwise interactions of the coordinates
#' mod2 <- semireg(log(grade) ~ s(zcoord, bs = "os") + s(xcoord, ycoord, bs = "ztps", k = 6) +
#'                   xcoord:ycoord + xcoord:zcoord + ycoord:zcoord, data = copper, prt = FALSE)
#' semipred(mod2, covariate = "zcoord", trans = exp, prt = FALSE)$plt
#' }
#' @aliases copperKnots
#' @name copper
#' @docType data
NULL

#' Growth of Indiana Children
#'
#' @description
#' Repeated height measurements of 216 adolescents, about every 6 months between the ages of 5 and 20, from a study
#' of the development of hypertension at the Indiana University School of Medicine, Indianapolis (Pratt et al.
#' 1989); the same data as \code{HRW::growthIndiana}. The data of M. P. Wand's examples of subject-specific
#' curves (Harezlak, Ruppert and Wand 2018, section 9.2), which \code{\link{semireg}} reproduces
#' (\code{tests/testthat/test-semireg-golden.R}).
#' @usage
#' data(growthIndiana)
#' @format
#' A data frame with 4123 rows and 5 columns:
#' \describe{
#' \item{idnum}{identification number of the adolescent.}
#' \item{height}{height (cm).}
#' \item{age}{age (years).}
#' \item{male}{1 for a boy, 0 for a girl.}
#' \item{black}{1 for a black adolescent, 0 otherwise.}
#' }
#' @source
#' Pratt, J. H., Jones, J. J., Miller, J. Z., Wagner, M. A. and Fineberg, N. S. (1989). Racial differences in
#' aldosterone excretion and plasma aldosterone concentrations in children. \emph{New England Journal of Medicine},
#' 321, 1152-1157.
#' @references
#' Harezlak, J., Ruppert, D. and Wand, M. P. (2018). \emph{Semiparametric Regression with R}. Springer.
#' @examples
#' data(growthIndiana)
#' # black boys: a mean curve, and each boy's own curve, a straight line (age | idnum)
#' # plus a spline deviation with one variance for all boys
#' gbm <- subset(growthIndiana, male == 1 & black == 1)
#' gbm$idnum <- factor(gbm$idnum)
#' fit <- semireg(height ~ s(age, bs = "os", k = 20, xt = list(range = c(5.5, 20))) +
#'                  (s(age, k = 10, xt = list(range = c(5.5, 20))) | idnum), data = gbm, prt = FALSE)
#' anova(fit)$smooth       # are the boys' curves straight lines?
#' semipred(fit, "age", "age", prt = FALSE)$plt
#' \donttest{
#' # all boys: a common curve and a black-white deviation curve (the contrast curve of
#' # Harezlak, Ruppert and Wand); at which ages do the races differ, and in growth rate?
#' gm <- subset(growthIndiana, male == 1)
#' gm$idnum <- factor(gm$idnum)
#' gm$race <- factor(ifelse(gm$black == 1, "black", "white"))
#' fit2 <- semireg(height ~ s(age, race, bs = "os", k = 20, xt = list(range = c(5.5, 20))) +
#'                   (s(age, k = 10, xt = list(range = c(5.5, 20))) | idnum), data = gm, prt = FALSE)
#' dif <- semipred(fit2, "race", "age", contr = c("black", "white"), band = "simultaneous", prt = FALSE)
#' attr(dif$pred_df, "significant")
#' semipred(fit2, "race", "age", deriv = 1, prt = FALSE)$plt     # growth rates
#' predictmeans(fit2, "race", covariate = c(age = 12))$means
#' }
#' @name growthIndiana
#' @docType data
NULL

#' Spinal Bone Mineral Density of Girls
#'
#' @description
#' Repeated measurements of spinal bone mineral density of 230 girls and young women aged 9 to 26 of four
#' ethnic groups, from a study of bone mineral acquisition (Bachrach et al. 1999); the data of M. P. Wand's
#' additive mixed model example (\code{gamm()}). \code{HRW::femSBMD} is a larger version (1003 measurements).
#' @usage
#' data(SBMDfem)
#' @format
#' A data frame with 547 rows and 7 columns:
#' \describe{
#' \item{idnum}{identification number of the subject.}
#' \item{spnbmd}{spinal bone mineral density (g/cm\eqn{^2}).}
#' \item{age}{age (years).}
#' \item{ethnicity}{factor: Asian, Black, Hispanic, White.}
#' \item{black, hispanic, white}{indicators (1 = yes) of the ethnic groups other than Asian.}
#' }
#' @source
#' Bachrach, L. K., Hastie, T., Wang, M.-C., Narasimhan, B. and Marcus, R. (1999). Bone mineral acquisition in
#' healthy Asian, Hispanic, Black, and Caucasian youth: a longitudinal study. \emph{Journal of Clinical
#' Endocrinology and Metabolism}, 84, 4702-4712.
#' @examples
#' data(SBMDfem)
#' SBMDfem$idnum <- factor(SBMDfem$idnum)
#' # a curve over age, ethnicity as a fixed effect, a random intercept per girl
#' fit <- semireg(spnbmd ~ s(age, bs = "os", k = 15) + ethnicity + (1 | idnum), data = SBMDfem, prt = FALSE)
#' anova(fit)       # ethnicity: F test with Kenward-Roger df; the smooth: likelihood ratio test
#' # the ethnic groups at age 15, compared in pairs
#' predictmeans(fit, "ethnicity", covariate = c(age = 15), pairwise = TRUE, adj = "tukey")$pairs
#' # growth of bone density: where is the curve rising?
#' attr(semipred(fit, "age", "age", deriv = 1, prt = FALSE)$pred_df, "significant")
#' @name SBMDfem
#' @docType data
NULL

#' Ozone Concentration and Meteorology in California
#'
#' @description
#' Daily ozone concentration with three meteorological variables (345 days), a version of the Los Angeles ozone
#' data of Breiman and Friedman (1985); the data of M. P. Wand's additive model examples (fitted by \code{lme()}
#' and by \code{gam()}).
#' @usage
#' data(ozoneCalif)
#' @format
#' A data frame with 345 rows and 4 columns:
#' \describe{
#' \item{ozone}{ozone concentration.}
#' \item{inversion.base.height}{inversion base height.}
#' \item{daggett.pressure.gradient}{pressure gradient (Los Angeles to Daggett).}
#' \item{inversion.base.temp}{inversion base temperature.}
#' }
#' @source
#' Breiman, L. and Friedman, J. H. (1985). Estimating optimal transformations for multiple regression and
#' correlation. \emph{Journal of the American Statistical Association}, 80, 580-598.
#' @examples
#' data(ozoneCalif)
#' # an additive model: three O'Sullivan splines
#' fit <- semireg(ozone ~ s(inversion.base.height, bs = "os", k = 15) +
#'                  s(daggett.pressure.gradient, bs = "os", k = 15) +
#'                  s(inversion.base.temp, bs = "os", k = 15), data = ozoneCalif, prt = FALSE)
#' anova(fit)$smooth     # is each curve more than a straight line?
#' # the effect of one covariate, the others at their means
#' semipred(fit, covariate = "inversion.base.temp", prt = FALSE)$plt
#' @name ozoneCalif
#' @docType data
NULL

#' Scallop Abundance off Long Island
#'
#' @description
#' Scallop catches at 148 locations of a 1990 survey cruise on the Atlantic continental shelf off Long Island,
#' New York (Ecker and Heltshe 1994; the data of \code{HRW::scallop}), with the knots of M. P. Wand's thin-plate
#' spline fit (\code{scallopKnots}) and the boundary of the survey area (\code{scallopBdry}).
#' @usage
#' data(scallop)
#' data(scallopKnots)
#' data(scallopBdry)
#' @format
#' \code{scallop}: a data frame with 148 rows and 3 columns:
#' \describe{
#' \item{latitude}{degrees latitude (north).}
#' \item{longitude}{degrees longitude (negative: west).}
#' \item{tot.catch}{total scallop catch at the location.}
#' }
#' \code{scallopKnots} (37 rows) and \code{scallopBdry} (77 rows, a polygon): data frames with columns
#' \code{longitude} and \code{latitude}.
#' @source
#' Ecker, M. D. and Heltshe, J. F. (1994). Geostatistical estimates of scallop abundance. In \emph{Case Studies in
#' Biometry}, eds. N. Lange et al., Wiley, New York, 107-124.
#' @references
#' Ruppert, D., Wand, M. P. and Carroll, R. J. (2003). \emph{Semiparametric Regression}. Cambridge University Press.
#' @examples
#' data(scallop)
#' data(scallopKnots)
#' data(scallopBdry)
#' # a thin-plate surface of asinh(catch) (low-rank kriging) with the given knots
#' fit <- semireg(asinh(tot.catch) ~ s(longitude, latitude, bs = "ztps",
#'                                    xt = list(knots = as.matrix(scallopKnots))),
#'                data = scallop, prt = FALSE)
#' anova(fit)$smooth
#' # the surface inside the survey area, on the scale of the catch
#' semipred(fit, covariate = c("longitude", "latitude"), boundary = scallopBdry, trans = sinh,
#'          prt = FALSE)$plt
#' @aliases scallopKnots scallopBdry
#' @name scallop
#' @docType data
NULL

#' Trade Union Membership and Wages
#'
#' @description
#' Hourly wage and trade union membership of 534 workers (the Current Population Survey of 1985, Berndt 1991),
#' the data of M. P. Wand's penalised spline logistic regression example.
#' @usage
#' data(wageTUnion)
#' @format
#' A data frame with 534 rows and 2 columns:
#' \describe{
#' \item{wage}{hourly wage (US dollars).}
#' \item{unionMember}{1 for a union member, 0 otherwise.}
#' }
#' @source
#' Berndt, E. R. (1991). \emph{The Practice of Econometrics: Classic and Contemporary}. Addison-Wesley.
#' @references
#' Ruppert, D., Wand, M. P. and Carroll, R. J. (2003). \emph{Semiparametric Regression}. Cambridge University Press.
#' @examples
#' data(wageTUnion)
#' # the probability of union membership as a smooth function of the wage
#' fit <- semireg(unionMember ~ s(wage, bs = "os", k = 15), family = binomial, data = wageTUnion, prt = FALSE)
#' anova(fit)$smooth
#' semipred(fit, "wage", "wage", prt = FALSE)$plt
#' @name wageTUnion
#' @docType data
NULL
