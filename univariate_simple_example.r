
# source the file
source("univariate.R")

# This data is the blood pressure example in winBUGS
bloodpressure<-list(sex = c(0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0,
1, 1, 1, 1),
age = c(59, 52, 37, 40, 67, 43, 61, 34, 51, 58, 54, 31, 49, 45, 66, 48, 41, 47, 53, 62,
60, 33, 44, 70, 56, 69, 35, 36, 68, 38),
bp = c(143, 132, 88, 98, 177, 102, 154, 83, 131, 150, 131, 69, 111, 114, 170, 117,
96, 116, 131, 158, 156, 75, 111, 184, 141, 182, 74, 87, 183, 89))

# turn data into create X and Y matrices
 X<- cbind(rep(1,30),bloodpressure$sex,bloodpressure$age)
 y <- bloodpressure$bp
# assign priors - we'll assume sigma2 is specified by the MEAN AND VARIANCE
 priorBeta = list(sigmaB=diag(c(1000,1000,1000)),B.init=c(-23,1,4))
 priorSigma2 = list(mu=1,variance=1,form="mean-var")   
 MCMCpars <- list(n.burn=1000,n.iter=20000,n.thin=1)
# do the sampling - use the known posteriors to sample directly
mod.1 <- univariate.all(y, X, priorB=priorBeta, priorSig=priorSigma2,MCMCpars = MCMCpars,makePlots=T) 
names(mod.1)


# Run a second model, using different priors
# assign priors - we'll assume sigma2 is specified by the MEAN AND VARIANCE
 priorSigma2 = list(n0=1,sig20=100,form="shape-scale")  
# do the sampling - use the known posteriors to sample directly
mod.2 <- univariate.vat(y, X, priorB=priorBeta, priorSig=priorSigma2,MCMCpars = MCMCpars,makePlots=T) 
names(mod.2)


##### Finally, show the same model using MCMCpack
# Advantages: blazing fast, does essentially the same computation
# Problems: 
library(MCMCpack)
example <- MCMCregress(bp ~ age+sex, data=bloodpressure,b0 = 0,
B0 = 1/100, c0 = 1, d0 = 1, marginal.likelihood = "Chib",burnin=10000,mcmc=100000)