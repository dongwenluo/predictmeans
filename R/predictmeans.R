predictmeans <- function (model, modelterm, data=NULL, pairwise=FALSE, atvar=NULL, adj="none", Df=NULL, lsd_bar=TRUE,
                          level=NULL, covariate=NULL, meandecr=NULL, letterCI=FALSE, trans = I, transOff=0,
						  responsen=NULL, count=FALSE, plotord=NULL, lineplot=TRUE, plottitle=NULL, plotxlab=NULL,
						  plotylab=NULL,
						  which.plots = c("mean", "pval", "back"),
						  mplot=NULL, barplot=NULL, pplot=NULL, bkplot=NULL, plot=TRUE, jitterv=0.2,
						  basesz=12, prtnum=TRUE, prtplt=TRUE, newwd=FALSE, permlist=NULL, ncore=3, ndecimal=4,
                          cov_reduce=mean, exclude=NULL, vcov_adjust=TRUE,
                          weights=c("equal", "proportional", "outer", "cells", "flat")) {
  vcov_adjust_default <- missing(vcov_adjust)

  ## Bias output towards fixed notation. See ?options for more info
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)

  if(any(missing(model), missing(modelterm))){
    stop("The arguments 'model', and 'modelterm' must be provided!")
  }

  weights <- match.arg(weights)

  # The calculation runs in stages that share this call's state (see
  # predictmeans-stages.R).
  st <- environment()
  pm_stage_setup(st)
  pm_stage_estimate(st)
  pm_stage_compare(st)
  pm_stage_mean_table(st)
  pm_stage_backtransform(st)
  if (plot) {
    pm_stage_plots(st)
  }
  pm_stage_assemble(st)
  get("pm_out", envir = st)
}
