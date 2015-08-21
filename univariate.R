# Code written by Eric Ward, 07.04.08
# email: eric.ward@noaa.gov
# 
# This is designed to do the univariate sampling for Bayesian linear regression
# Two versions are included, one utilizing multivariate sampling, one doing variable at a time.
# Output from this program can be used as an mcmc() object, and fed directly into CODA for
# diagnostics, by calling codamenu().  
#
# Sources:        Unless otherwise specified, all eqns and page refs refer to Gelman
# 1. Gelman (2004), Bayesian Data Analysis 2nd ed
# 2. Carlin and Louis () Bayes and Empirical Bayes methods for data analysis
# 3. Chib (1995) Marginal Likelihood from the Gibbs Output, JASA


library(coda)     # needed if mcmc objects are returned
library(mvtnorm)  # needed to call dmvnorm
log2pi <- log(2*pi)
#library(MASS)     # needed to call rmvnorm
# Functions for random number generation, not currently used
#rinvgamma<-function(shape, scale = 1) return(1/rgamma(1, shape, scale))
#rinvchi2<-function(shape, scale = 1) return(shape*scale/rchisq(1,shape))
#rinvchi2.2 <-function(shape, scale = 1) return(shape*scale/rchisq(1,shape))
# Function for estimating prior density
dinvgamma.log <- function(x, shape, scale = 1) 
{
    log.density <- shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - (scale/x)
    return(log.density)
}
dinvchi2.log <- function(x, shape, scale = 1) return(dinvgamma.log(x, shape/2, shape*scale/2))
dinvchi2 <- function(x, shape, scale = 1) return(exp(dinvgamma.log(x, shape/2, shape*scale/2)))      
# Function to compute MV normal random numbers, mvrnorm in MASS takes 137% longer, because it uses eigen 
rmvnorm.chol <- function(norm.01, u, Sigma) {
   # norm.01 is a vector of Normal(0,1) random numbers
   return(c(t(u) + (norm.01)%*%chol(Sigma)))
}


###############################################################################
#  Function for doing sampling posteriors for univariate normal linear regression, 
#  with conjugate priors.  The difference between univariate.exact and univariate.gibbs
#  is that univariate.exact draws all B parameters from a MVN distribution.
###############################################################################
univariate.all <- function(y, X, priorB = list(B0=NA, sigmaB=NA,B.init=NA), priorSig = list(form="shape-scale",n0=NA,sig20=NA,mu=NA,variance=NA,sig.init=NA), MCMCpars=list(n.burn=1000, n.iter=10000, n.thin=1),makePlots=FALSE) {
  # y is the data, can be passed in as a vector or matrix
  # X is a matrix of covariates, including column of 1s for intercept
  # priorB is an optional list containing: 
  #    B0, vector of prior means (defaults to 0)
  #    sigmaB, vector / matrix of prior covariances (defaults to diagonal)
  #    B.init, vector of initial values (defaults to 0)
  # priorSig is an optional list containing: 
  #    form, either "shape-scale" (default) or "mean-var"
  #    if the prior shape-scale (default): 
  #    n0: how many data points is this prior worth? (defaults to 1)
  #    sig20, prior scale parameter (defaults to 1) 
  #    if the prior mean-var:
  #    mu, prior mean parameter
  #    variance, prior var parameter
  #    sig.init, initial value (defaults to MLE)
  # MCMCpars is an optional list containing:
  #    n.burn, period of burn in (true length = burn * thin)
  #    n.iter, number of samples to save after the burn
  #    n.thin, thinning interval
  # makePlots is T/F, whether to plot the log(BF) and MCMC object, default = F

  # the "form" argument of the priorSig list allows the user to optionally specify the prior for sigma2
  # in terms of the mean and the variance
  if(priorSig$form=="mean-var") {
    # then mean and var need to be translated into shape / scale
    priorSig$n0 = (2 * priorSig$mu^2)/priorSig$variance + 4        # shape parameter
    priorSig$sig20 = (priorSig$mu)*(priorSig$n0-2)/priorSig$n0     # scale parameter
  }

  n = length(y)    # number of observations
  k = dim(X)[2]    # number of estimated pars, not including sig2
  # The next 4 variables are just constants to be used later
  n.min.k = n-k
  n0plusnk = priorSig$n0 + (n + k)
  nplusk = n+k
  n0sig20 = priorSig$n0+priorSig$sig20
  #if(is.na(sigmaY)==TRUE) sigmaY = diag(n)  else sigmaY = priorSig$sigmaY # if sigmaY is not passed in, use identity matrix
  Qy = diag(n)   # this is scaling cov matrix, sigmaY = Qy * sigma2
  # Make sure that the priors for B have been specified correctly
  if(is.null(priorB$B0) == TRUE || length(priorB$B0) != k || any(is.na(priorB$B0)==TRUE)) priorB$B0 = rep(0,k)         # default prior has mean 0
  if(is.null(priorB$B.init) == TRUE || length(priorB$B.init) != k || any(is.na(priorB$B.init)==TRUE)) priorB$B.init = rep(0,k)
  B = priorB$B.init # Assign B to the initial guesses               
  
  # Make sure that the priors for sigma2 have been specified correctly
  #if(priorSig$n0 <= 4) stop("Error: increase shape parameter for sigma2 prior")
  #if(priorSig$sig20 <= 0) stop("Error: increase scale parameter for sigma2 prior")
  if(is.null(priorSig$sig.init) == TRUE || is.na(priorSig$sig.init) == TRUE) {
      Bhat = chol2inv(chol(t(X)%*%X))%*%t(X)%*%y
      priorSig$sig.init = sqrt(t(y-X%*%Bhat)%*%(y-X%*%Bhat)/(n-k))
  }
  sigma2 = priorSig$sig.init
  # Create X.star and Y.star, eq 14.25 Gelman p 384
  X.star = rbind(X, diag(k))
  Y.star = as.matrix(c(y,priorB$B0))
  # These matrices are zeros to fill in the components of Sigma.star
  zeros.UR = matrix(0,n,k)
  zeros.LL = matrix(0,k,n)                                       
  
  if(is.null(priorB$sigmaB) == TRUE || dim(as.matrix(priorB$sigmaB))[1] != k || any(is.na(priorB$sigmaB)==TRUE)==TRUE) priorB$sigmaB = rep(1,k)
  if(is.matrix(priorB$sigmaB)==TRUE) diag.priorB = (priorB$sigmaB) else diag.priorB = diag(priorB$sigmaB)

  #################################################################
  # THIS IS THE BURN-IN PHASE
  #################################################################
  invSigma.star = diag(1/diag(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))
  rndchi2 = rchisq(n=MCMCpars$n.burn*MCMCpars$n.thin, df=n0plusnk)  # take advantage of vectors, draw all the Chi2 draws we'll need
  for(i in 1:(MCMCpars$n.burn*MCMCpars$n.thin)) {
     # Y | B, sig2 ~ norm(XB, sig2*I)
     # sig2 | y ~ InvX2(n0+n+k, (n0*sig20+(n+k)*s2)/(n0+n+k)  # n+k is length of Y.star
     y.min.XB = Y.star - X.star%*%B                           # residuals
     s2 = t(y.min.XB)%*%invSigma.star%*%y.min.XB/(n.min.k)          # eqn 14.18
     # draw x ~ chi2 (v) RV, v*s2/x ~ InvChi2(v, s2)          # Inv Chi2, Appendix, p580
     sigma2 = (n0sig20+(nplusk)*s2)/rndchi2[i]                # p384 
     # Multivariate distribution of B | sigma2, y
     # MVN (Bhat, Vb)  Bhat = inv(t(X)*X) *t(X)*y, Vb = inv(t(X)*X)
     # We can build the star matrices (14.25), and then plug them in as known covariances (14.11-14.12)
     # this is for all cases, EW commented out the next 2 lines for speed in the independent case
     # invSigma.star = chol2inv(chol(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))   
     # avoid inverting the matrix because we can assume the data & pars are independent
     invSigma.star = diag(1/diag(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))
     # Use cholesky decomp to avoid calculating inverse
     Vb = chol2inv(chol(t(X.star)%*%invSigma.star%*%X.star))  # eqn 14.3 or 14.11
     Bhat = Vb%*%t(X.star)%*%invSigma.star%*%Y.star           # 14.12 
     B = rmvnorm.chol(rnorm(k,0,1),Bhat,Vb*c(sigma2))         # 14.3
  }  
  
  #################################################################
  # THIS IS THE SAMPLING PHASE
  #################################################################
  c = 1
  rndchi2 = rchisq(n=MCMCpars$n.iter*MCMCpars$n.thin, df=n0plusnk)  # take advantage of vectors, draw all the Chi2 draws we'll need
  B.post = matrix(0,MCMCpars$n.iter,k)
  sigma2.post = 0
  for(i in 1:MCMCpars$n.iter) {    # Nested for loops avoid the if or modulus statements - makes speed 30-50% faster
     for(j in 1:MCMCpars$n.thin) {
        # Y | B, sig2 ~ norm(XB, sig2*I)
        # sig2 | y ~ InvX2(n0+n+k, (n0*sig20+(n+k)*s2)/(n0+n+k)  # n+k is length of Y.star
        y.min.XB = Y.star - X.star%*%B                           # residuals
        s2 = t(y.min.XB)%*%invSigma.star%*%y.min.XB/(n.min.k)          # eqn 14.18
        # draw x ~ chi2 (v) RV, v*s2/x ~ InvChi2(v, s2)          # Inv Chi2, Appendix, p580
        sigma2 = (n0sig20+(nplusk)*s2)/rndchi2[c]                # p384
        c = c+1 
        # Multivariate distribution of B | sigma2, y
        # MVN (Bhat, Vb)  Bhat = inv(t(X)*X) *t(X)*y, Vb = inv(t(X)*X)
        # We can build the star matrices (14.25), and then plug them in as known covariances (14.11-14.12)
        # this is for all cases, EW commented out the next 2 lines for speed in the independent case
        # invSigma.star = chol2inv(chol(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))   
        # avoid inverting the matrix because we can assume the data & pars are independent
        invSigma.star = diag(1/diag(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))
        # Use cholesky decomp to avoid calculating inverse
        Vb = chol2inv(chol(t(X.star)%*%invSigma.star%*%X.star))  # eqn 14.3 or 14.11
        Bhat = Vb%*%t(X.star)%*%invSigma.star%*%Y.star           # 14.12 
        B = rmvnorm.chol(rnorm(k,0,1),Bhat,Vb*c(sigma2))         # 14.3
     }
     # Store the parameter draws
     B.post[i,] = B
     sigma2.post[i] = sigma2
  }    

  #################################################################
  # THIS IS THE BAYES FACTOR CALCULATION / PARAMETER SUMMARY PHASE
  #################################################################  
  B.mode = 0        
  for(i in 1:k) {
     d = density(B.post[,i])
     B.mode[i] = d$x[which(d$y==max(d$y))]   # Calculate the mode
  }
  d = density(sqrt(sigma2.post))
  sigma.mode = d$x[which(d$y==max(d$y))]     # Calculate the mode
  sigma2.mode <- sigma.mode^2
  # calculate the log likelihood at the posterior mode
  logLike = sum(dnorm(x = y, mean = X%*%B.mode, sd = sigma.mode,log=T))  # calculate the log likelihood
  # calculate log prior for sigma2, alternatively could estimate mode for tau and use gamma  
  logPrior = dinvchi2.log(x = sigma.mode^2,shape = priorSig$n0, scale = priorSig$sig20)  # Appendix A
  logPrior = logPrior + dmvnorm(x=B.mode, mean = priorB$B0, sigma=diag.priorB, log=TRUE)
  
  # log(posterior) is broken down into two points,
  # Posterior(B*,sigma*|Y) = Posterior(B*|sigma*,Y)*Posterior(sigma*|Y)
  # See Chib, section 2.1.2 or maybe a better explanation in Carlin & Louis, p 209
  # Posterior(B*|sigma*,Y) is directly available:
  sigma2 <- sigma2.mode
  invSigma.star = diag(1/diag(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))
  Vb = chol2inv(chol(t(X.star)%*%invSigma.star%*%X.star))  # eqn 14.3 or 14.11
  Bhat = Vb%*%t(X.star)%*%invSigma.star%*%Y.star           # 14.12 
  logPost.B <- dmvnorm(x=B.mode, mean = Bhat, sigma=Vb*c(sigma2), log=TRUE)
  logPost.sigma <- 0
  for(i in 1:MCMCpars$n.iter) {  # Rao-Blackwell approximation, Carlin & Louis, p 209 eqn 6.1.6
     # Calculate Posterior(sigma*|Y) as average over P(sigma*|B[i],Y) 
     y.min.XB = Y.star - X.star%*%B.post[i,]         # residuals
     invSigma.star = diag(1/diag(rbind(cbind(Qy * c(sigma2.mode),zeros.UR),cbind(zeros.LL,diag.priorB))))
     s2 = t(y.min.XB)%*%invSigma.star%*%y.min.XB/(n.min.k)          # eqn 14.18
     logPost.sigma[i] <- dinvchi2(x = sigma2.mode, shape = n0plusnk, scale = s2)
  }
  logPost.sigma <- log(mean(logPost.sigma))
  
  # marginal likelihood - essentially log [p(model | data) ]
  # Chib refers to it as BMI, basic marginal likelihood identity
  BFML = logLike + logPrior - logPost.B - logPost.sigma  # Chib, section 2.1.2  
  
  # return data frame of saved values
  post = as.data.frame(cbind(B.post,sqrt(sigma2.post)))
  for(i in 1:k) names(post)[i] = paste("B",i-1,sep="")
  names(post)[k+1] = "sigma"
  print(summary(mcmc(post)))    # summarize output to console
  print(paste("Log(Bayes Factor): ",BFML,sep=""),quote=FALSE)
  ############  THIS IS TO MAKE OPTIONAL PLOTS
  if(makePlots == TRUE) {
     #par(ask=T) 
     plot(mcmc(post))         
     
     par(mfrow=c(2,2)) # make 2 x 2 panel plots of priors & posteriors
     #num.plots = ceil((k+1)/4)
     Bpriors <- matrix(0,20000,k)
     for(i in 1:20000) Bpriors[i,] = rmvnorm.chol(rnorm(k,0,1),priorB$B0, diag.priorB)
     for(i in 1:k) {
        # draw a sample from the prior variance to plot
        post = B.post[,i]
        prior = Bpriors[,i]
        #min.x = min(c(sigma2.post)); max.x = max(c(sigma2.post))
        min.x = min(post,prior); max.x = max(post,prior)
        d.prior = density(prior); d.prior$y = d.prior$y/sum(d.prior$y)
        d.post = density(post); d.post$y = d.post$y/sum(d.post$y) 
        #min.x = min(d.post$x[which(d.post$y > median(d.post$y))])
        #max.x = max(d.post$x[which(d.post$y > median(d.post$y))])
        #d.like = 
        min.y = 0.05*max(d.prior$y,d.post$y);  max.y = max(1.10*max(d.prior$y,d.post$y),0.01+max(d.prior$y,d.post$y))
        plot(d.prior$x,d.prior$y,main="",ylab="Density",xlab=paste("B",i-1,sep=""),xlim=c(min.x,max.x),ylim=c(min.y,max.y),type="l",col="blue",lwd=2)
        #lines(d.like$x,d.like$y,col="red",lwd=2)
        lines(d.post$x,d.post$y,lwd=2,col="purple")
        text(d.prior$x[which(d.prior$y==max(d.prior$y))], 0.004+max(d.prior$y),expression(paste("P(",theta,")")),cex=0.7)
        #text(d.like$x[which(d.like$y==max(d.like$y))], 0.006+max(d.like$y),expression(paste("L(Y|",theta,")")))
        text(d.post$x[which(d.post$y==max(d.post$y))], 0.004+max(d.post$y),expression(paste("P(",theta,"|Y)")),cex=0.7)     
     }
     
     # draw a sample from the prior variance to plot
     prior = log((priorSig$n0*priorSig$sig20)/rchisq(100000,df=priorSig$n0)) 
     post = sigma2.post
     #min.x = min(c(sigma2.post)); max.x = max(c(sigma2.post))
     min.x = min(post,prior); max.x = max(post,prior)
     d.prior = density(prior,from=0); d.prior$y = d.prior$y/sum(d.prior$y)
     d.post = density(post,from=0); d.post$y = d.post$y/sum(d.post$y) 
     #min.x = min(d.post$x[which(d.post$y > median(d.post$y))])
     #max.x = max(d.post$x[which(d.post$y > median(d.post$y))])
     #d.like = 
     min.y = 0.05*max(d.prior$y,d.post$y);  max.y = max(1.10*max(d.prior$y,d.post$y),0.01+max(d.prior$y,d.post$y))
     plot(d.prior$x,d.prior$y,main="",ylab="Density",xlab=expression(sigma[2]),xlim=c(min.x,max.x),ylim=c(min.y,max.y),type="l",col="blue",lwd=2)
     #lines(d.like$x,d.like$y,col="red",lwd=2)
     lines(d.post$x,d.post$y,lwd=2,col="purple")
     text(d.prior$x[which(d.prior$y==max(d.prior$y))], 0.004+max(d.prior$y),expression(paste("P(",theta,")")),cex=0.7)
     #text(d.like$x[which(d.like$y==max(d.like$y))], 0.006+max(d.like$y),expression(paste("L(Y|",theta,")")))
     text(d.post$x[which(d.post$y==max(d.post$y))], 0.004+max(d.post$y),expression(paste("P(",theta,"|Y)")),cex=0.7)  
  }
  list(saved = post, B.mode = B.mode,sigma.mode=sigma.mode,LogBayesFac=BFML)
}



###############################################################################
#  Function for doing variable at a time (vat) gibbs sampling for univariate normal
#  linear regression, with conjugate priors
###############################################################################

univariate.vat <- function(y, X, priorB = list(B0=NA, sigmaB=NA,B.init=NA), priorSig = list(form="shape-scale",n0=NA,sig20=NA,mu=NA,variance=NA,sig.init=NA), MCMCpars=list(n.burn=1000,n.iter=10000,n.thin=1),makePlots=FALSE) {
  # y is the data, can be passed in as a vector or matrix
  # X is a matrix of covariates, including column of 1s for intercept
  # priorB is an optional list containing: 
  #    B0, vector of prior means (defaults to 0)
  #    sigmaB, vector / matrix of prior covariances (defaults to diagonal)
  #    B.init, vector of initial values (defaults to 0)
  # priorSig is an optional list containing: 
  #    form, either "shape-scale" (default) or "mean-var"
  #    if the prior shape-scale (default): 
  #    n0, prior shape parameter: how many data points is this prior worth? (defaults to 1)
  #    sig20, prior scale parameter (defaults to 1) 
  #    if the prior mean-var:
  #    mu, prior mean parameter
  #    variance, prior var parameter
  #    sig.init, initial value (defaults to MLE)
  # MCMCpars is an optional list containing:
  #    n.burn, period of burn in (true length = burn * thin)
  #    n.iter, number of samples to save after the burn
  #    n.thin, thinning interval
  # makePlots is T/F, whether to plot the log(BF) and MCMC object, default = F

  # the "form" argument of the priorSig list allows the user to optionally specify the prior for sigma2
  # in terms of the mean and the variance
  if(priorSig$form=="mean-var") {
    # then mean and var need to be translated into shape / scale
    priorSig$n0 = (2 * priorSig$mu^2)/priorSig$variance + 4        # shape parameter
    priorSig$sig20 = (priorSig$mu)*(priorSig$n0-2)/priorSig$n0     # scale parameter
  }

  n = length(y)    # number of observations
  k = dim(X)[2]    # number of estimated pars, not including sig2
  # The next 4 variables are just constants to be used later
  n.min.k = n-k
  n0plusnk = priorSig$n0 + (n + k)
  nplusk = n+k
  n0sig20 = priorSig$n0+priorSig$sig20
  #if(is.na(sigmaY)==TRUE) sigmaY = diag(n)  else sigmaY = priorSig$sigmaY # if sigmaY is not passed in, use identity matrix
  Qy = diag(n)   # this is scaling cov matrix, sigmaY = Qy * sigma2
  # Make sure that the priors for B have been specified correctly
  if(is.null(priorB$B0) == TRUE || length(priorB$B0) != k || any(is.na(priorB$B0)==TRUE)) priorB$B0 = rep(0,k)         # default prior has mean 0
  if(is.null(priorB$B.init) == TRUE || length(priorB$B.init) != k || any(is.na(priorB$B.init)==TRUE)) priorB$B.init = rep(0,k)
  B = priorB$B.init # Assign B to the initial guesses               
  
  # Make sure that the priors for sigma2 have been specified correctly
  #if(priorSig$n0 <= 4) stop("Error: increase shape parameter for sigma2 prior")
  #if(priorSig$sig20 <= 0) stop("Error: increase scale parameter for sigma2 prior")
  if(is.null(priorSig$sig.init) == TRUE || is.na(priorSig$sig.init) == TRUE) {
      Bhat = chol2inv(chol(t(X)%*%X))%*%t(X)%*%y
      priorSig$sig.init = sqrt(t(y-X%*%Bhat)%*%(y-X%*%Bhat)/(n-k))
  }
  sigma2 = priorSig$sig.init
  # Create X.star and Y.star, eq 14.25 Gelman p 384
  X.star = rbind(X, diag(k))
  Y.star = as.matrix(c(y,priorB$B0))
  # These matrices are zeros to fill in the components of Sigma.star
  zeros.UR = matrix(0,n,k)
  zeros.LL = matrix(0,k,n)                                       
  
  if(is.null(priorB$sigmaB) == TRUE || dim(as.matrix(priorB$sigmaB))[1] != k || any(is.na(priorB$sigmaB)==TRUE)==TRUE) priorB$sigmaB = rep(1,k)
  if(is.matrix(priorB$sigmaB)==TRUE) diag.priorB = (priorB$sigmaB) else diag.priorB = diag(priorB$sigmaB)

  #################################################################
  # THIS IS THE BURN-IN PHASE
  #################################################################
  invSigma.star = diag(1/diag(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))
  rndchi2 = rchisq(n=MCMCpars$n.burn*MCMCpars$n.thin, df=n0plusnk)  # take advantage of vectors, draw all the Chi2 draws we'll need
  for(i in 1:(MCMCpars$n.burn*MCMCpars$n.thin)) {
     # Y | B, sig2 ~ norm(XB, sig2*I)
     # sig2 | y ~ InvX2(n0+n+k, (n0*sig20+(n+k)*s2)/(n0+n+k)  # n+k is length of Y.star
     y.min.XB = Y.star - X.star%*%B                           # residuals
     s2 = t(y.min.XB)%*%invSigma.star%*%y.min.XB/(n.min.k)          # eqn 14.18
     # draw x ~ chi2 (v) RV, v*s2/x ~ InvChi2(v, s2)          # Inv Chi2, Appendix, p580
     sigma2 = (n0sig20+(nplusk)*s2)/rndchi2[i]                # p384 
     # Multivariate distribution of B | sigma2, y
     # MVN (Bhat, Vb)  Bhat = inv(t(X)*X) *t(X)*y, Vb = inv(t(X)*X)
     # We can build the star matrices (14.25), and then plug them in as known covariances (14.11-14.12)
     # this is for all cases, EW commented out the next 2 lines for speed in the independent case
     # invSigma.star = chol2inv(chol(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))   
     # avoid inverting the matrix because we can assume the data & pars are independent
     invSigma.star = diag(1/diag(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))
     # Use cholesky decomp to avoid calculating inverse
     Vb = chol2inv(chol(t(X.star)%*%invSigma.star%*%X.star))  # eqn 14.3 or 14.11
     Bhat = Vb%*%t(X.star)%*%invSigma.star%*%Y.star           # 14.12 
     # Do variable at a time gibbs sampling
     Vb.inv = chol2inv(chol(Vb*c(sigma2)))  # calculate inverese of cov matrix
     ind.vars = (1/diag(Vb.inv))        # this is used in the var argument of A.2, and in sq brackets on line above
     C = diag(k)-diag(ind.vars)%*%Vb.inv # calculate conditional coefficients, eqn A.2 in Appendix
     for(j in 1:k) {
        # sample the conditional univariate distribution, A.2 in Appendix
        # Mean can be computed several ways: 2nd method preferred for speed
        #   mu[j] = Bhat[j] + Vb.sigma2[j,-j]%*%chol2inv(chol(Vb.sigma2[-j,-j]))%*%as.matrix((B[-j]-Bhat[-j]))
        #   mu[j] = Bhat[j] + C[j,]%*%as.matrix(B-Bhat)
        # Variance can be computed several ways: 2nd method preferred for speed
        #   var[j] = Vb.sigma2[j,j] - Vb.sigma2[j,-j]%*%chol2inv(chol(Vb.sigma2[-j,-j]))%*%Vb.sigma2[-j,j]
        #   var[j] = ind.vars[j]
        B[j] = rnorm(1, Bhat[j] + C[j,]%*%as.matrix(B-Bhat), sqrt(ind.vars[j]))
     }     
  }  
  
  #################################################################
  # THIS IS THE SAMPLING PHASE
  #################################################################
  c = 1
  rndchi2 = rchisq(n=MCMCpars$n.iter*MCMCpars$n.thin, df=n0plusnk)  # take advantage of vectors, draw all the Chi2 draws we'll need
  B.post = matrix(0,MCMCpars$n.iter,k)
  sigma2.post = 0
  for(i in 1:MCMCpars$n.iter) {    # Nested for loops avoid the if or modulus statements - makes speed 30-50% faster
     for(j in 1:MCMCpars$n.thin) {
        # Y | B, sig2 ~ norm(XB, sig2*I)
        # sig2 | y ~ InvX2(n0+n+k, (n0*sig20+(n+k)*s2)/(n0+n+k)  # n+k is length of Y.star
        y.min.XB = Y.star - X.star%*%B                           # residuals
        s2 = t(y.min.XB)%*%invSigma.star%*%y.min.XB/(n.min.k)          # eqn 14.18
        # draw x ~ chi2 (v) RV, v*s2/x ~ InvChi2(v, s2)          # Inv Chi2, Appendix, p580
        sigma2 = (n0sig20+(nplusk)*s2)/rndchi2[c]                # p384
        c = c+1 
        # Multivariate distribution of B | sigma2, y
        # MVN (Bhat, Vb)  Bhat = inv(t(X)*X) *t(X)*y, Vb = inv(t(X)*X)
        # We can build the star matrices (14.25), and then plug them in as known covariances (14.11-14.12)
        # this is for all cases, EW commented out the next 2 lines for speed in the independent case
        # invSigma.star = chol2inv(chol(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))   
        # avoid inverting the matrix because we can assume the data & pars are independent
        invSigma.star = diag(1/diag(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))
        # Use cholesky decomp to avoid calculating inverse
        Vb = chol2inv(chol(t(X.star)%*%invSigma.star%*%X.star))  # eqn 14.3 or 14.11
        Bhat = Vb%*%t(X.star)%*%invSigma.star%*%Y.star           # 14.12 
        # Do variable at a time gibbs sampling
        Vb.inv = chol2inv(chol(Vb*c(sigma2)))  # calculate inverese of cov matrix
        ind.vars = (1/diag(Vb.inv))        # this is used in the var argument of A.2, and in sq brackets on line above
        C = diag(k)-diag(ind.vars)%*%Vb.inv # calculate conditional coefficients, eqn A.2 in Appendix
        for(j in 1:k) {
           # sample the conditional univariate distribution, A.2 in Appendix
           # Mean can be computed several ways: 2nd method preferred for speed
           #   mu[j] = Bhat[j] + Vb.sigma2[j,-j]%*%chol2inv(chol(Vb.sigma2[-j,-j]))%*%as.matrix((B[-j]-Bhat[-j]))
           #   mu[j] = Bhat[j] + C[j,]%*%as.matrix(B-Bhat)
           # Variance can be computed several ways: 2nd method preferred for speed
           #   var[j] = Vb.sigma2[j,j] - Vb.sigma2[j,-j]%*%chol2inv(chol(Vb.sigma2[-j,-j]))%*%Vb.sigma2[-j,j]
           #   var[j] = ind.vars[j]
           B[j] = rnorm(1, Bhat[j] + C[j,]%*%as.matrix(B-Bhat), sqrt(ind.vars[j]))
        }     
     }
     # Store the parameter draws
     B.post[i,] = B
     sigma2.post[i] = sigma2
  }    
  
  #################################################################
  # THIS IS THE BAYES FACTOR CALCULATION / PARAMETER SUMMARY PHASE
  #################################################################  
  B.mode = 0        
  for(i in 1:k) {
     d = density(B.post[,i])
     B.mode[i] = d$x[which(d$y==max(d$y))]   # Calculate the mode
  }
  d = density(sqrt(sigma2.post))
  sigma.mode = d$x[which(d$y==max(d$y))]     # Calculate the mode
  sigma2.mode <- sigma.mode^2
  # calculate the log likelihood at the posterior mode
  logLike = sum(dnorm(x = y, mean = X%*%B.mode, sd = sigma.mode,log=T))  # calculate the log likelihood
  # calculate log prior for sigma2, alternatively could estimate mode for tau and use gamma  
  logPrior = dinvchi2.log(x = sigma.mode^2,shape = priorSig$n0, scale = priorSig$sig20)  # Appendix A
  logPrior = logPrior + dmvnorm(x=B.mode, mean = priorB$B0, sigma=diag.priorB, log=TRUE)
  
  # log(posterior) is broken down into two points,
  # Posterior(B*,sigma*|Y) = Posterior(B*|sigma*,Y)*Posterior(sigma*|Y)
  # See Chib, section 2.1.2 or maybe a better explanation in Carlin & Louis, p 209
  # Posterior(B*|sigma*,Y) is directly available:
  sigma2 <- sigma2.mode
  invSigma.star = diag(1/diag(rbind(cbind(Qy * c(sigma2),zeros.UR),cbind(zeros.LL,diag.priorB))))
  Vb = chol2inv(chol(t(X.star)%*%invSigma.star%*%X.star))  # eqn 14.3 or 14.11
  Bhat = Vb%*%t(X.star)%*%invSigma.star%*%Y.star           # 14.12 
  logPost.B <- dmvnorm(x=B.mode, mean = Bhat, sigma=Vb*c(sigma2), log=TRUE)
  logPost.sigma <- 0
  for(i in 1:MCMCpars$n.iter) {  # Rao-Blackwell approximation, Carlin & Louis, p 209 eqn 6.1.6
     # Calculate Posterior(sigma*|Y) as average over P(sigma*|B[i],Y) 
     y.min.XB = Y.star - X.star%*%B.post[i,]         # residuals
     invSigma.star = diag(1/diag(rbind(cbind(Qy * c(sigma2.mode),zeros.UR),cbind(zeros.LL,diag.priorB))))
     s2 = t(y.min.XB)%*%invSigma.star%*%y.min.XB/(n.min.k)          # eqn 14.18
     logPost.sigma[i] <- dinvchi2(x = sigma2.mode, shape = n0plusnk, scale = s2)
  }
  logPost.sigma <- log(mean(logPost.sigma))
  
  # marginal likelihood - essentially log [p(model | data) ]
  # Chib refers to it as BMI, basic marginal likelihood identity
  BFML = logLike + logPrior - logPost.B - logPost.sigma  # Chib, section 2.1.2  
  
  # return data frame of saved values
  post = as.data.frame(cbind(B.post,sqrt(sigma2.post)))
  for(i in 1:k) names(post)[i] = paste("B",i-1,sep="")
  names(post)[k+1] = "sigma"
  print(summary(mcmc(post)))    # summarize output to console
  print(paste("Log(Bayes Factor): ",BFML,sep=""),quote=FALSE)
  ############  THIS IS TO MAKE OPTIONAL PLOTS
  if(makePlots == TRUE) {
     #par(ask=T) 
     plot(mcmc(post))         
     
     par(mfrow=c(2,2)) # make 2 x 2 panel plots of priors & posteriors
     #num.plots = ceil((k+1)/4)
     Bpriors <- matrix(0,20000,k)
     for(i in 1:20000) Bpriors[i,] = rmvnorm.chol(rnorm(k,0,1),priorB$B0, diag.priorB)
     for(i in 1:k) {
        # draw a sample from the prior variance to plot
        post = B.post[,i]
        prior = Bpriors[,i]
        #min.x = min(c(sigma2.post)); max.x = max(c(sigma2.post))
        min.x = min(post,prior); max.x = max(post,prior)
        d.prior = density(prior); d.prior$y = d.prior$y/sum(d.prior$y)
        d.post = density(post); d.post$y = d.post$y/sum(d.post$y) 
        #min.x = min(d.post$x[which(d.post$y > median(d.post$y))])
        #max.x = max(d.post$x[which(d.post$y > median(d.post$y))])
        #d.like = 
        min.y = 0.05*max(d.prior$y,d.post$y);  max.y = max(1.10*max(d.prior$y,d.post$y),0.01+max(d.prior$y,d.post$y))
        plot(d.prior$x,d.prior$y,main="",ylab="Density",xlab=paste("B",i-1,sep=""),xlim=c(min.x,max.x),ylim=c(min.y,max.y),type="l",col="blue",lwd=2)
        #lines(d.like$x,d.like$y,col="red",lwd=2)
        lines(d.post$x,d.post$y,lwd=2,col="purple")
        text(d.prior$x[which(d.prior$y==max(d.prior$y))], 0.004+max(d.prior$y),expression(paste("P(",theta,")")),cex=0.7)
        #text(d.like$x[which(d.like$y==max(d.like$y))], 0.006+max(d.like$y),expression(paste("L(Y|",theta,")")))
        text(d.post$x[which(d.post$y==max(d.post$y))], 0.004+max(d.post$y),expression(paste("P(",theta,"|Y)")),cex=0.7)     
     }
     
     # draw a sample from the prior variance to plot
     prior = log((priorSig$n0*priorSig$sig20)/rchisq(100000,df=priorSig$n0)) 
     post = sigma2.post
     #min.x = min(c(sigma2.post)); max.x = max(c(sigma2.post))
     min.x = min(post,prior); max.x = max(post,prior)
     d.prior = density(prior,from=0); d.prior$y = d.prior$y/sum(d.prior$y)
     d.post = density(post,from=0); d.post$y = d.post$y/sum(d.post$y) 
     #min.x = min(d.post$x[which(d.post$y > median(d.post$y))])
     #max.x = max(d.post$x[which(d.post$y > median(d.post$y))])
     #d.like = 
     min.y = 0.05*max(d.prior$y,d.post$y);  max.y = max(1.10*max(d.prior$y,d.post$y),0.01+max(d.prior$y,d.post$y))
     plot(d.prior$x,d.prior$y,main="",ylab="Density",xlab=expression(sigma[2]),xlim=c(min.x,max.x),ylim=c(min.y,max.y),type="l",col="blue",lwd=2)
     #lines(d.like$x,d.like$y,col="red",lwd=2)
     lines(d.post$x,d.post$y,lwd=2,col="purple")
     text(d.prior$x[which(d.prior$y==max(d.prior$y))], 0.004+max(d.prior$y),expression(paste("P(",theta,")")),cex=0.7)
     #text(d.like$x[which(d.like$y==max(d.like$y))], 0.006+max(d.like$y),expression(paste("L(Y|",theta,")")))
     text(d.post$x[which(d.post$y==max(d.post$y))], 0.004+max(d.post$y),expression(paste("P(",theta,"|Y)")),cex=0.7)  
  }
  list(saved = post, B.mode = B.mode,sigma.mode=sigma.mode,LogBayesFac=BFML)
}