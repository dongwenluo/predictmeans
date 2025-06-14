#' \code{ATP} containing data
#'
#' \code{ATP} containing data from an experiment to study the effects of
#' preserving liquids on the enzyme content of dog hearts. There were 23 hearts
#' and two treatment factors, A and B, each at two levels. Measurements were
#' made of ATP as a percentage of total enzyme in the heart, at one and two
#' hourly intervals during a twelve hour period following initial preservation.
#'
#'
#' @name ATP
#' @docType data
#' @format \code{ATP} is a 230 row data frame with the following columns
#' \describe{ \item{heart}{dog heart id.} \item{time}{time in hour for
#' \code{ATP} measurement.} \item{A}{treatment with two levels.}
#' \item{B}{treatment with two levels.} \item{ATP}{percentage of total enzyme
#' in the heart.} }
NULL





#' \code{Clinical} data
#'
#' \code{Clinical} data is from a multicentre randomized clinical trial
#' (Beitler & Landis 1985, Biometrics).
#'
#'
#' @name Clinical
#' @docType data
#' @format \code{Clinical} is a 16 row data frame with the following columns
#' \describe{ \item{Clinic}{8 centres id.} \item{Treatment}{2 skin treatments
#' (control or active drug).} \item{Favorable}{number that produced a
#' favourable response.} \item{Total}{number of patients in each treatment
#' group.} }
NULL





#' \code{Drug} data
#'
#' The data is for the comparison of the effectiveness of three analgesic drugs
#' to a standard drug, morphine (Finney, Probit analysis, 3rd Edition 1971,
#' p.103). 14 groups of mice were tested for response to the drugs at a range
#' of doses.
#'
#'
#' @name Drug
#' @docType data
#' @format \code{Drug} is a 14 row data frame with the following columns
#' \describe{ \item{Drug}{type of drug.} \item{Dose}{dose volumn.}
#' \item{N}{total number of mice in each group.} \item{R}{number responding
#' mice in each group.} \item{log10Dose}{log10 transformed dose volumn.} }
NULL





#' Predicted Means for Linear and Semiparametric Models
#'
#' This package provides functions to diagnose and make inferences from various
#' linear models, such as those obtained from 'aov', 'lm', 'glm', 'gls', 'lme',
#' 'lmer', 'glmmTMB' and 'semireg'. Inferences include predicted means and
#' standard errors, contrasts, multiple comparisons, permutation tests,
#' adjusted R-square and graphs.
#'
#' \tabular{ll}{ Package: \tab predictmeans\cr Type: \tab Package\cr Version:
#' \tab 1.1.2\cr Date: \tab 2025-05-08\cr License: \tab GPL (>= 2)\cr }
#'
#' @name predictmeans-package
#' @author Dongwen Luo, James Curran, Simon Urbanek, Siva Ganesh and John
#' Koolaard
#'
#' Maintainer: Dongwen Luo <dongwen.luo@@agresearch.co.nz>
#' @references Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R.
#' (2004), \emph{Prediction in linear mixed models}, Australian and New Zealand
#' Journal of Statistics, 46(3), 325-347.
#' @keywords internal
"_PACKAGE"
