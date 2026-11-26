## Help pages of the data sets (roxygen2 sources).

#' ATP Content of Preserved Dog Hearts
#'
#' @description
#' Data from an experiment to study the effects 
#' of preserving liquids on the enzyme content of dog hearts. There were 23 hearts 
#' and two treatment factors, A and B, each at two levels. Measurements were made 
#' of ATP as a percentage of total enzyme in the heart, at one and two hourly
#' intervals during a twelve hour period following initial preservation.
#' @usage
#' data(ATP)
#' @format
#' A data frame with 230 rows (23 hearts, 10 times) and 5 columns:
#' \describe{
#' \item{heart}{factor, 23 dog hearts.}
#' \item{time}{factor, time of measurement (hours after preservation; 10 levels).}
#' \item{A}{factor, treatment A (2 levels).}
#' \item{B}{factor, treatment B (2 levels).}
#' \item{ATP}{ATP as a percentage of total enzyme in the heart (\%).}
#' }
#' @examples
#' data(ATP)
#' fm <- lmer(ATP ~ A*B*time + (1|heart), data = ATP)
#' predictmeans(fm, "A:B", pairwise = TRUE)$means
#' @name ATP
#' @docType data
NULL

#' Body Length of Beetles (Simulated)
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
#' R2_model(fm)
#' @name Body_length
#' @docType data
NULL

#' Multicentre Clinical Trial of a Skin Treatment
#'
#' @description
#' Favourable responses in a multicentre randomized clinical trial of a topical drug against a control (Beitler
#'  & Landis 1985, Biometrics).
#' @usage
#' data(Clinical)
#' @format
#' A data frame with 16 rows and 4 columns:
#' \describe{
#' \item{Clinic}{factor, 8 centres.}
#' \item{Treatment}{factor, \code{control} or \code{drug}.}
#' \item{Favorable}{number of patients with a favourable response.}
#' \item{Total}{number of patients in the clinic and treatment group.}
#' }
#' @source
#' Beitler, P. J. and Landis, J. R. (1985). A mixed-effects model for categorical data. \emph{Biometrics}, 41, 991-1000.
#' @examples
#' data(Clinical)
#' gm <- glmer(cbind(Favorable, Total - Favorable) ~ Treatment + (1|Clinic),
#'             family = binomial, data = Clinical)
#' predictmeans(gm, "Treatment")$means   # probability scale
#' @name Clinical
#' @docType data
NULL

#' Potency of Analgesic Drugs in Mice
#'
#' @description
#' Data for the comparison of the effectiveness of three analgesic drugs 
#' to a standard drug, morphine (Finney, Probit analysis, 3rd Edition 1971, p.103). 14 groups 
#' of mice were tested for response to the drugs at a range of doses.
#' @usage
#' data(Drug)
#' @format
#' A data frame with 14 rows and 5 columns:
#' \describe{
#' \item{Drug}{factor, 4 drugs (Morphine, Amidone, ...).}
#' \item{Dose}{dose (units as in Finney, 1971).}
#' \item{N}{number of mice in the group.}
#' \item{R}{number of mice responding.}
#' \item{log10Dose}{\code{log10(Dose)}.}
#' }
#' @source
#' Finney, D. J. (1971). \emph{Probit Analysis}, 3rd edition. Cambridge University Press, p. 103.
#' @examples
#' data(Drug)
#' gm <- glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial(link = "probit"), data = Drug)
#' predictmeans(gm, "Drug", covariate = c(log10Dose = 0.5))$means
#' @name Drug
#' @docType data
NULL

#' Egg Counts of Beetles (Simulated)
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
#' @examples
#' data(Poisson_ind)
#' gm <- glmer(Egg ~ Treatment + Habitat + (1|Population/Container) + (1|individual:Container),
#'             family = poisson, data = Poisson_ind)
#' R2_model(gm)
#' R2_model(gm, over_disp = TRUE)   # without the overdispersion term
#' @name Poisson_ind
#' @docType data
NULL

#' Copper Grade at Sampling Points
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
#' \item{grade}{copper grade (positive; modelled on the log scale below).}
#' }
#' \code{copperKnots}: the 20 knots (columns \code{xcoord}, \code{ycoord}) of M. P. Wand's thin-plate spline
#' fit of these data (a geoadditive model fitted by \code{lme()}), which \code{\link{semireg}} reproduces.
#' @examples
#' data(copper)
#' data(copperKnots)
#' # M. P. Wand's geoadditive model: a thin-plate surface over the plane with the given knots,
#' # and an O'Sullivan spline in depth (their linear terms are added to the fixed effects)
#' mod <- semireg(grade ~ s(xcoord, ycoord, bs = "ztps", xt = list(knots = as.matrix(copperKnots))) +
#'                  s(zcoord, bs = "os", xt = list(knots = seq(80, 120, by = 5))), data = copper)
#' anova(mod)   # fixed effects (Kenward-Roger) and a test of each smooth term
#' semipred(mod, covariate = c("xcoord", "ycoord"))$plt
#' \donttest{
#' # on the log scale, with the pairwise interactions of the coordinates
#' mod2 <- semireg(log(grade) ~ s(zcoord, bs = "os") + s(xcoord, ycoord, bs = "ztps", k = 6) +
#'                   xcoord:ycoord + xcoord:zcoord + ycoord:zcoord, data = copper)
#' semipred(mod2, covariate = "zcoord", trans = exp)$plt
#' }
#' @aliases copperKnots
#' @name copper
#' @docType data
NULL

#' Heights of Trees: 37 Seedlots at 6 Sites
#'
#' @description
#' Mean heights of 37 seedlots (26 Acacia, Allocasuarina and Casuarina species) grown at six sites in Thailand,
#' the data of Table 5.2 of Williams, Harwood and Matheson (2024), used there for a regression of each seedlot on
#' the site mean (Finlay-Wilkinson regression). 30 site-seedlot combinations have no height.
#' @usage
#' data(DataExam5.2)
#' @format
#' A data frame with 222 rows and 4 columns:
#' \describe{
#' \item{site}{factor, 6 sites.}
#' \item{seedlot}{factor, 37 seedlot numbers.}
#' \item{ht}{height (cm); NA when the seedlot was not measured at the site.}
#' \item{sitemean}{the mean of \code{ht} at the site (one value per site).}
#' }
#' @details
#' \code{sitemean} takes one value per site, so in \code{ht ~ site + sitemean * seedlot} it is aliased with
#' \code{site} (its coefficient is NA); see the section "Covariates aliased with a factor" of
#' \code{\link{predictmeans}}. On the four sites and 27 seedlots with complete data (and \code{sitemean}
#' recomputed for them) \code{predictmeans()} reproduces the means and slopes of the book's Table 5.7
#' (\code{tests/testthat/test-level-covariate.R}).
#' @source
#' Williams, E. R., Harwood, C. E. and Matheson, A. C. (2024). \emph{Experimental Design and Analysis for Tree
#' Improvement} (3rd ed.). CSIRO Publishing, Table 5.2. \doi{10.1071/9781486317110}
#' @examples
#' data(DataExam5.2)
#' fit <- lm(ht ~ site + sitemean * seedlot, data = DataExam5.2)
#' head(predictmeans(fit, "seedlot", plot = FALSE)$mean_table)
#' # each seedlot averaged over the sites where it was measured
#' head(predictmeans(fit, "seedlot", plot = FALSE, weights = "cells")$mean_table)
#' @name DataExam5.2
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
#'                  (s(age, k = 10, xt = list(range = c(5.5, 20))) | idnum), data = gbm)
#' anova(fit)$smooth       # are the boys' curves straight lines?
#' semipred(fit, "age", "age")$plt
#' \dontrun{
#' # (about 30 s: 2257 measurements, 130 boys' curves)
#' # all boys: a common curve and a black-white deviation curve (the contrast curve of
#' # Harezlak, Ruppert and Wand); at which ages do the races differ, and in growth rate?
#' gm <- subset(growthIndiana, male == 1)
#' gm$idnum <- factor(gm$idnum)
#' gm$race <- factor(ifelse(gm$black == 1, "black", "white"))
#' fit2 <- semireg(height ~ s(age, race, bs = "os", k = 20, xt = list(range = c(5.5, 20))) +
#'                   (s(age, k = 10, xt = list(range = c(5.5, 20))) | idnum), data = gm)
#' dif <- semipred(fit2, "race", "age", contr = c("black", "white"), band = "simultaneous")
#' attr(dif$pred_df, "significant")
#' semipred(fit2, "race", "age", deriv = 1)$plt     # growth rates
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
#' fit <- semireg(spnbmd ~ s(age, bs = "os", k = 15) + ethnicity + (1 | idnum), data = SBMDfem)
#' anova(fit)       # ethnicity: F test with Kenward-Roger df; the smooth: likelihood ratio test
#' # the ethnic groups at age 15, compared in pairs
#' predictmeans(fit, "ethnicity", covariate = c(age = 15), pairwise = TRUE, adj = "tukey")$pairs
#' # growth of bone density: where is the curve rising?
#' attr(semipred(fit, "age", "age", deriv = 1)$pred_df, "significant")
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
#'                  s(inversion.base.temp, bs = "os", k = 15), data = ozoneCalif)
#' anova(fit)$smooth     # is each curve more than a straight line?
#' # the effect of one covariate, the others at their means
#' semipred(fit, covariate = "inversion.base.temp")$plt
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
#'                data = scallop)
#' anova(fit)$smooth
#' # the surface inside the survey area, on the scale of the catch
#' semipred(fit, covariate = c("longitude", "latitude"), boundary = scallopBdry, trans = sinh)$plt
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
#' fit <- semireg(unionMember ~ s(wage, bs = "os", k = 15), family = binomial, data = wageTUnion)
#' anova(fit)$smooth
#' semipred(fit, "wage", "wage")$plt
#' @name wageTUnion
#' @docType data
NULL

#' California Test Score Data
#'
#' @description
#' Test scores, school characteristics and student demographics of 420 school districts in California, 1998-1999
#' (kindergarten to grade 6 or 8), as used by Stock and Watson (2007) and in package \pkg{Ecdat}. The data of the
#' semiparametric example with several smooth terms and two random intercepts in the overview article
#' (\code{vignettes/articles/overview.Rmd}).
#' @usage
#' data(Caschool)
#' @format
#' A data frame with 420 rows and 17 columns:
#' \describe{
#' \item{distcod}{district code.}
#' \item{county}{county (a factor with 45 levels).}
#' \item{district}{district (a factor).}
#' \item{grspan}{grade span of the district: \code{"KK-06"} or \code{"KK-08"}.}
#' \item{enrltot}{total enrollment.}
#' \item{teachers}{number of teachers.}
#' \item{calwpct}{percentage of students qualifying for CalWORKs (income assistance).}
#' \item{mealpct}{percentage of students qualifying for a reduced-price lunch.}
#' \item{computer}{number of computers.}
#' \item{testscr}{average test score (the mean of \code{readscr} and \code{mathscr}).}
#' \item{compstu}{computers per student.}
#' \item{expnstu}{expenditure per student.}
#' \item{str}{student-teacher ratio.}
#' \item{avginc}{district average income (in USD 1000).}
#' \item{elpct}{percentage of English learners.}
#' \item{readscr}{average reading score.}
#' \item{mathscr}{average math score.}
#' }
#' @source
#' California Department of Education, data of 1998-1999; package \pkg{Ecdat} (Croissant and Graves), data set
#' \code{Caschool}.
#' @references
#' Stock, J. H. and Watson, M. W. (2007). \emph{Introduction to Econometrics}, 2nd edition. Addison-Wesley,
#' chapters 4-7.
#' @examples
#' data(Caschool)
#' Caschool$log.avginc <- log(Caschool$avginc)
#' \donttest{
#' # math scores against three covariates, each a smooth curve, with random
#' # intercepts for counties and grade spans
#' fit <- semireg_tmb(mathscr ~ s(calwpct, bs = "os", k = 12) + s(log.avginc, bs = "os", k = 5) +
#'                      s(compstu, bs = "os", k = 12) + (1 | county) + (1 | grspan), data = Caschool)
#' semipred(fit, covariate = "log.avginc")$plt
#' }
#' @name Caschool
#' @docType data
NULL
