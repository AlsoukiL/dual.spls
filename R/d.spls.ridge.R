#' Dual Sparse Partial Least Squares (Dual-SPLS) regression for the norm \eqn{\Omega(w)=\lambda \|w\|_1 +lambda2 \|Xw\|_2 + \|w\|_2}
#' @description
#' The function \code{d.spls.ridge} performs dimentional reduction combined to variable selection using the
#' Dual-SPLS algorithm with the norm \eqn{\Omega(w)=\lambda \|w\|_1 +lambda2 \|Xw\|_2 + \|w\|_2}.
#' @usage d.spls.ridge(X,y,ncp,pctnu,nu2,verbose=FALSE)
#' @param X a numeric matrix of predictors values. Each row represents an observation and each column a predictor variable.
#' @param y a numeric vector or a one column matrix of responses. It represents the response variable for each converstation.
#' @param ncp a positive integer. \code{ncp} is the number of Dual-SPLS components.
#' @param pctnu a positive real value, in \code{[0,1]}. \code{pctnu} is the desired
#' proportion of variables to shrink to zero for each component.
#' @param nu2 a positive real value. \code{nu2} is a sparsity parameter.
#' @param verbose a boolean value indicating whether or not to diplay the iterations steps.
#' @details
#' This procedure computes latent sparse components that are used in a regression model.
#' The optimization problem of the Dual-SPLS regression for the norm \eqn{\Omega(w)=\lambda \|w\|_1 +lambda2 \|Xw\|_2 + \|w\|_2} is
#' similar to the Dual norm defintion of the norm \eqn{\Omega(w)}. Indeed, we are searching for \code{w} that goes
#' with
#' \deqn{\Omega^*(z)=max_w(z^Tw) \text{ s.t. } \Omega(w)=1}
#' The solution of this problem is
#' \deqn{\dfrac{w_j}{\|w\|_2}=\dfrac{1}{\mu} \delta_{N_{2,j}} (|z_{N_{2},\nu_2}|-\nu_1)_+}
#' Where
#' \itemize{
#' \item \eqn{\delta_{N_{2,j}}} is the vector of signs of the product between \eqn{(\nu_2 N_2^TN_2+I_p)^{-1}} and z.
#' \item \eqn{\mu} is a parameter that guarentees the constraint of \eqn{\Omega(w)=1}
#' \item \eqn{\nu_1} is the shrinkage parameter.
#' \item \eqn{z_{N_{2},\nu_2}} is the product between \eqn{(\nu_2 N_2^TN_2+I_p)^{-1}} and z.
#'
#' }
#'
#' The parameter \eqn{\nu_2} is computed adaptively according to the proportion of zero variables desired \code{pctnu}.
#' @return A \code{list} of the following attributes
#' \itemize{
#' \item  Xmean the mean vector of the predictors matrix \code{X}.
#' \item  scores a matrix of \code{ncp} columns representing the Dual-SPLS components and the same number of rows
#' as \code{X} representing the observations in the new component basis computed by the compression step
#' of the Dual-SPLS.
#' \item  loadings the loadings matrix.
#' \item  Bhat the regression coefficients matrix for each component.
#' \item  intercept the intercept value for each component.
#' \item  fitted.values the matrix of predicted values of \code{y}
#' \item  residuals the matrix of residuals corresponding to the difference between the response and the fitted values.
#' \item  lambda1 the vector of the sparse hyper-parameter used to fit the model at each iteration.
#' \item  zerovar the vector of the number of variables shrinked to zero per component.
#' }
#' @author François Wahl Louna Alsouki
#' @seealso `browseVignettes("dual.spls")`
#'
#' @examples
#' ### load dual.spls library
#' library(dual.spls)
#' ### parameters
#' n <- 200
#' p <- 100
#' nondes <- 150
#' sigmaondes <- 0.01
#' data.benchmark=BCHMK(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
#'
#' X <- data.benchmark$X
#' y <- data.benchmark$y
#'
#' #fitting the model
#' mod.dspls <- d.spls.ridge(X=X,y=y,ncp=10,pctnu=0.9,nu2=0.05,verbose=TRUE)
#'
#' str(mod.dspls)
#'
#' ### plotting the observed values VS predicted values
#' plot(y,mod.dspls$fitted.values[,6], xlab="Observed values", ylab="Predicted values",
#'  main="Observed VS Predicted for 6 components")
#' points(-1000:1000,-1000:1000,type='l')
#'
#' ### plotting the regression coefficients
#' par(mfrow=c(3,1))
#' i=6
#' nz=mod.dspls$zerovar[i]
#' plot(1:dim(X)[2],mod.dspls$Bhat[,i],type='l',
#'     main=paste(" Dual-SPLS (ridge), ncp =", i, " #0coef =", nz, "/", dim(X)[2]),
#'     ylab='',xlab='' )
#' inonz=which(mod.dspls$Bhat[,i]!=0)
#' points(inonz,mod.dspls$Bhat[inonz,i],col='red',pch=19,cex=0.5)
#' legend("topright", legend ="non null values", bty = "n", cex = 0.8, col = "red",pch=19)
#' @export


d.spls.ridge<- function(X,y,ncp,pctnu,nu2,verbose=FALSE)
{

  ###################################
  # Dimensions
  ###################################
  n=length(y) #Number of observations
  p=dim(X)[2] #Number of variables

  ###################################
  # Centering Data
  ###################################
  Xm = apply(X, 2, mean) #Mean of X
  Xc=X - rep(1,n) %*% t(Xm) #Centering predictor matrix

  ym=mean(y) #Mean of y
  yc=y-ym #Centering response vector

  ###################################
  # Initialisation
  ###################################
  WW=matrix(0,p,ncp) #Initialising W, the matrix of loadings
  TT=matrix(0,n,ncp) #Initialising T, the matrix of scores
  Bhat=matrix(0,p,ncp) #Initialising the matrix of coefficients
  YY=matrix(0,n,ncp) #Initialising the matrix of coefficients
  RES=matrix(0,n,ncp) #Initialising the matrix of coefficients
  intercept=rep(0,ncp)
  zerovar=rep(0,ncp)
  listelambda1=rep(0,ncp)

  ###################################
  # Dual-SPLS
  ###################################

  Xdef=Xc #Initialising X for Deflation Step
  #Each step ic in -for loop- determine the icth column of each W, T and Bhat
  for (ic in 1:ncp)
  {

    # Z=X^Ty
    z=t(Xdef)%*%yc #For cov(t(X)y,w)=0, w must be colinear to t(X)y ==> Z=t(X)y
    z=as.vector(z)

    # Computing z12

    #Computing (nu2 N2T N2 + I)^(-1)
    N2=Xdef
    N2TN2=t(Xdef)%*%Xdef
    inv=solve(nu2*N2TN2+diag(p))


    #Computing delta
    delta=sign(z)

    ZN=inv%*%z

    #Optimizing nu
    zs=sort(abs(ZN))
    zsp=(1:p)/p
    iz=which.min(abs(zsp-pctnu))
    ###########
    nu1=zs[iz] #
    # ###########

    # finding lambda, mu, given nu
    z12=sapply(ZN,function(u) sign(u)*max(abs(u)-nu1,0))

    # Computing mu

    #########
    mu=norm2(z12) #
    ##############
    lambda1=nu1/mu #
    ##############

    # calculating w,t at the optimum
    w=(mu/(nu2*norm2(N2%*%z12)^2+nu1*norm1(z12) + mu^2))*z12

    #Finding T
    t=Xdef%*%w
    t=t/norm2(t)

    WW[,ic]=w
    TT[,ic]=t

    #Deflation
    Xdef=Xdef-t%*%t(t)%*%Xdef

    #Coefficient vectors
    R=t(TT[,1:ic,drop=FALSE])%*%Xc%*%WW[,1:ic,drop=FALSE]
    R[row(R)>col(R)]<-0 # inserted for numerical stability

    L=backsolve(R,diag(ic))
    Bhat[,ic]=WW[,1:ic,drop=FALSE]%*%(L%*%(t(TT[,1:ic,drop=FALSE])%*%yc))

    listelambda1[ic]=lambda1
    intercept[ic] = ym - Xm %*% Bhat[,ic]
    zerovar[ic]=sum(Bhat[,ic]==0)

    #Predictions
    YY[,ic]=X %*% Bhat[,ic] + intercept[ic]
    RES[,ic]=y-YY[,ic]

    # results iteration
    if (verbose){
      cat('Dual PLS ic=',ic,'lambda1=',lambda1,
          'mu=',mu,'nu2=',nu2,
          'nbzeros=',zerovar[ic], '\n')
    }

  }

  return(list(Xmean=Xm,scores=TT,loadings=WW,Bhat=Bhat,intercept=intercept,
              fitted.values=YY,residuals=RES,
              lambda1=listelambda1,zerovar=zerovar))

}
