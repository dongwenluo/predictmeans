########## R script: SBMDfemGCV ##########

# For analysis of the spinal bone mineral
# density data using the function gamm() in
# in the package mgcv.

# Last changed: 11 APR 2011 by M.P. Wand

# Load required packages:

library(lattice) ; library(mgcv) 

# Read in and extract data to be analysed:

SBMDfem <- read.table("SBMDfem.txt",header=TRUE)

# Obtain plot of raw data:

SBMDfemVis <- xyplot(spnbmd~age|factor(ethnicity),group=idnum,
                as.table=TRUE,data=SBMDfem,xlab="age (years)",
                ylab=expression(paste("spinal bone mineral
                density (g/c",m^2,")")),
                panel=function(x,y,subscripts,groups)
                {  
                   panel.grid() 
                   panel.superpose(x,y,subscripts,groups,
                                   type="b",pch=16,lwd=2)
                })
print(SBMDfemVis)
wait()

# Do additive mixed model using gamm():

fit <- gamm(spnbmd~s(age)+black+hispanic+white,
            random=list(idnum=~1),data=SBMDfem)


# Plot the fitted age curve:

plot(fit$gam,shade=TRUE,shade.col="palegreen")

# Obtain a residual plot:

plot(fit$lme)

# Obtain numerical summaries:

print(summary(fit$gam)) ; wait()
print(intervals(fit$lme)) ; wait()

# Do likelihood ratio test for ethnicity effect:

fitReduced <- gamm(spnbmd~s(age),random=list(idnum=~1),data=SBMDfem)
logLRstat <- 2*(fit$lme$logLik - fitReduced$lme$logLik)
pValue <- 1 - pchisq(logLRstat,df=3)
cat("Approx. p-value for likelihood ratio test for ethnicity effect is:\n")
print(pValue)

############ End of SBMDfemGCV ##########
