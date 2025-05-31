#' Permutation ANOVA for \code{lmer} Model
#'
#' This function provides permutation ANOVA for \code{lmer} model.
#'
#'
#' @param model Model object returned by \code{lmer}.
#' @param nperm Number of permutation, the default value is 999.
#' @param ncore Number of core for parallel computing, the default value is 3.
#' @param type The type of ANOVA table requested (using SAS terminology) with
#' Type I being the familiar sequential ANOVA table.
#' @param ...  Use to setup option: seed -- Specify a random number generator
#' seed, for reproducible results.
#' @return Permutation ANOVA table.
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'
#' \dontrun{
#'  library(predictmeans)
#'  Oats$nitro <- factor(Oats$nitro)
#'  fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'
#'  # Permutation Test for model terms
#'  permanova.lmer(fm)
#'  permanova.lmer(fm, type=2)
#'
#'  # Compare to F test
#'  fm0 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#'  anova(fm0)
#' }
#' @importFrom lme4 refit
#' @importFrom stats anova
#' @importFrom utils as.roman
#' @export
permanova.lmer <- function(model, nperm = 999, ncore=3L, type = c("I", "II", "III",  "1", "2", "3"), ...) {

  if (!inherits(model, "lmerMod")) {
    stop("The model must be a lmer object!")
  }
  if (inherits(try(refit(model, getME(model, "y")), TRUE), "try-error")) {
    stop("Please remove missing values in your model!")
  }
  type <- type[1L]
  if (!is.character(type)) {
    type <- as.character(type)
  }
  type <- match.arg(type)
  if (type %in% c("I", "II")) {
    type <- as.character(as.integer(as.roman(type)))
  }

  aTable <- anova(model, type=type)
  Perm.p <- numeric(nrow(aTable))
  termlabel1 <- termlabel2 <- row.names(aTable)
  names(Perm.p) <- termlabel1

  termlabel0 <- attr(terms(model), "term.labels")
  model.0 <- update(model, as.formula(paste(".~. -", paste(termlabel0, collapse="-"))))

  for (vars in termlabel1) {
    if (type == "3") {
      varsn <- unlist(strsplit(vars, "\\:"))
      for (i in varsn) termlabel2 <- termlabel2[grep(i, termlabel2)]
      termlabel <- paste(termlabel2, collapse="-")
      model.b <- update( model, as.formula(paste(".~. -", termlabel)))
      Perm.p[vars] <- permlmer(model.b, model, nperm, ncore, plot=FALSE, ...)$`Perm-p`[2]
      termlabel2 <- termlabel1
    } else if (type=="2") {
      varsn <- unlist(strsplit(vars, "\\:"))
      model.1 <- update(model.0, as.formula(paste(".~.+", paste(termlabel0[attr(terms(model), "order")%in%(1:length(varsn))], collapse="+"))))
      model.2 <- update(model.1, as.formula(paste(".~. -", vars)))
      Perm.p[vars] <- permlmer(model.2, model.1, nperm, ncore, plot=FALSE, ...)$`Perm-p`[2]
    } else if (type=="1") {
      model.1 <- update(model.0, as.formula(paste(".~. +", vars)))
      Perm.p[vars] <- permlmer(model.0, model.1, nperm, ncore, plot=FALSE, ...)$`Perm-p`[2]
      model.0 <- model.1
    } else {
      NULL
    }
  }
  aTable$Perm.p <- Perm.p
  return(aTable)
}
