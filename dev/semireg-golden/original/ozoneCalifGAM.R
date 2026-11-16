########## R script: ozoneCalifGAM.Rs ##########

# For fitting a generalised additive model to
# ozone data from California, USA. In this case
# Generalised Cross Validation (GCV) is used
# for smoothing parameter selection.

# Last changed: 05 AUG 2010 by M.P. Wand.

# Load required package:

library(mgcv) 
source("wait.r")

# Read in, extract data and display:

ozoneCalif <- read.table("ozoneCalif.txt",header=T)
pairs(ozoneCalif)
wait()

# Fit generalised additive model and plot estimated functions:

fit <- gam(ozone~s(inversion.base.height)+s(daggett.pressure.gradient)
                 +s(inversion.base.temp),data=ozoneCalif)
par(mfrow=c(1,3))
plot(fit,se=TRUE,shade=TRUE,shade.col="cyan")

########## End of ozoneCalifGAM ##########


