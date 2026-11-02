CookD <- function (model, group=NULL, plot=TRUE, idn=3, newwd=FALSE) {

  stopifnot("CookD doesn't support this model!"={
    any(inherits(model, "gls"), inherits(model, "lme"), inherits(model, "lmerMod"))
  })
  if (inherits(model, "gls") || inherits(model, "lme")) {
    mdf <- as.data.frame(pm_get_data(model))
    model <- pm_refit(model, data=mdf, method="ML", check=TRUE)
  }
  if (inherits(model, "lmerMod")) {
    model <- pm_refit(model, REML=FALSE, check=TRUE)
    mdf <- model.frame(model)
  }
  mp <- mymodelparm(model, df = 0)
  beta0  <- mp$coef
  vcovb <- mp$vcov
  vb.inv <- solve(vcovb)

  if (is.null(group) || group%in%c("NULL", "")) {
    rn <- rownames(mdf)
    LOOmp <- lapply(rn, function(x)  mymodelparm(pm_refit(model, data=mdf[rn!=x, ]), df = 0))
  } else {
    rn <- unique(mdf[, group])
    LOOmp <- lapply(rn, function(x)  {
  	rind <- mdf[, group]!=x
  	mymodelparm(pm_refit(model, data=mdf[rind, ]), df = 0)
    })
  }

  LOObeta <- sapply(LOOmp, function(x) x$coef)          # Matrix with beta(-i)
  rK <- t(LOObeta-beta0)
  CookD <- diag(rK %*% tcrossprod(vb.inv, rK)/length(beta0))
  names(CookD) <- rn
  if (plot) {
    if (newwd) {
      dev.new()
    }
    outD <- CookD >= sort(CookD, decreasing =TRUE)[idn]                        # Outlying Di's
    labid <- names(CookD)
    plot(CookD,  xlab="Obs. number", col="blue", ylim=c(0, max(CookD)+0.005),
         main="Cook's Distance", ylab = "Cook's distance", type = "h")
    text((1:length(CookD))[outD], CookD[outD], labid[outD], pos=3)          # Annotation
    points((1:length(CookD))[outD], CookD[outD], pch=16, col="blue")
  }
  return(invisible(CookD))
}
