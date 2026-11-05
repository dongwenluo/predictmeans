# Help page copper (roxygen2 source; man/copper.Rd is generated).
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
