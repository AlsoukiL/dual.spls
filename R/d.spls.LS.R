#' Dual Sparse Partial Least Squares (Dual-SPLS) regression for the norm \eqn{\Omega(w)=\lambda \|N_1w\|_1 + \|Xw\|_2}
#' @description
#' The function \code{d.spls.lasso} performs dimensional reduction as in PLS methodology combined to variable selection via the
#' Dual-SPLS algorithm with the norm \eqn{\Omega(w)=\lambda \|N_1w\|_1 + \|Xw\|_2}.
#' @usage d.spls.LS(X,y,ncp,ppnu,verbose=FALSE)
#' @param X a numeric matrix of predictors values. Each row represents an observation and each column a predictor variable.
#' @param y a numeric vector or a one column matrix of responses. It represents the response variable for each observation.
#' @param ncp a positive integer. \code{ncp} is the number of Dual-SPLS components.
#' @param ppnu a positive real value, in \code{[0,1]}. \code{ppnu} is the desired
#' proportion of variables to shrink to zero for each component (see Dual-SPLS methodology).
#' @param verbose a boolean value indicating whether or not to diplay the iterations steps.
#' @details
#' This procedure computes latent sparse components that are used in a regression model.
#' The optimization problem of the Dual-SPLS regression for the norm \eqn{\Omega(w)=\lambda \|N_1w\|_1 + \|Xw\|_2}
#' comes from the Dual norm defintion of \eqn{\Omega(w)}. Indeed, we are searching for \code{w} that goes
#' with
#' \deqn{\Omega^*(z)=max_w(z^Tw) \text{ s.t. } \Omega(w)=1}
#' Noting that \eqn{\lambda} is the initial shrinkage parameter that imposes sparsity, the Dual-SPLS does not rely
#' on the value of \eqn{\lambda} but instead proceeds adaptively by choosing the propotion of zeros that the user
#' would like to impose in the coefficients. Which leads to the compuation of a secondary shrinkage parameter \eqn{\nu}.
#'
#' The solution of this problem is
#' \deqn{\dfrac{w_j}{\|Xw\|_2}=\dfrac{1}{\mu} sign(\hat{\beta}_{LS_j}) (|\hat{\beta}_{LS_j}|-\nu)_+}
#' Where
#' \itemize{
#' \item \eqn{\hat{\beta}_{LS_j}} is equal to \eqn{(X^TX)^{-1}z}, the classical Least Squares regression coefficients.
#' \item \eqn{\mu} is a parameter that guarentees the constraint of \eqn{\Omega(w)=1}
#' \item \eqn{\nu} is the shrinkage parameter.
#' }
#' @return A \code{list} of the following attributes
#' \item{Xmean}{the mean vector of the predictors matrix \code{X}.}
#' \item{scores}{the matrix of dimension \eqn{n x ncp} where \code{n} is the number of observations.The \code{scores} represents
#' the observations in the new component basis computed by the compression step
#' of the Dual-SPLS.}
#' \item{loadings}{the matrix of dimension \eqn{p x ncp} that represents the Dual-SPLS components.}
#' \item{Bhat}{the matrix of dimension \eqn{p x ncp} that regroups the regression coefficients for each component.}
#' \item{intercept}{the vector of intercept values for each component.}
#' \item{fitted.values}{the matrix of dimension \eqn{n x ncp} that represents the predicted values of \code{y}}
#' \item{residuals}{the matrix of dimension \eqn{n x ncp} that represents the residuals corresponding
#'  to the difference between the responses and the fitted values.}
#' \item{zerovar}{the vector of length \eqn{ncp} represnting the number of variables shrinked to zero per component.}
#' @author Louna Alsouki François Wahl
#' @seealso `browseVignettes("dual.spls")`
#'
#' @examples
#' ### load dual.spls library
#' library(dual.spls)
#' ### parameters
#' n <- 100
#' p <- 50
#' nondes <- 20
#' sigmaondes <- 0.5
#' data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
#'
#' X <- data$X
#' y <- data$y
#'
#' #fitting the model
#' mod.dspls <- d.spls.LS(X=X,y=y,ncp=10,ppnu=0.9,verbose=TRUE)
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
#'
#'i=6
#' nz=mod.dspls$zerovar[i]
#' plot(1:dim(X)[2],mod.dspls$Bhat[,i],type='l',
#'     main=paste(" Dual-SPLS (LS), ncp =", i, " #0coef =", nz, "/", dim(X)[2]),
#'     ylab='',xlab='' )
#' inonz=which(mod.dspls$Bhat[,i]!=0)
#' points(inonz,mod.dspls$Bhat[inonz,i],col='red',pch=19,cex=0.5)
#' legend("topright", legend ="non null values", bty = "n", cex = 0.8, col = "red",pch=19)
#' @export


d.spls.LS<- function(X,y,ncp,ppnu,verbose=FALSE)
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
  YY=matrix(0,n,ncp) #Initialising the matrix of estimated responses
  RES=matrix(0,n,ncp) #Initialising the matrix of residues
  intercept=rep(0,ncp)
  zerovar=rep(0,ncp)

  ###################################
  # Dual-SPLS
  ###################################

  #Each step ic in -for loop- determine the icth column of each W, T and Bhat
  Xi=Xc #Initialising X for Deflation Step
  for (ic in 1:ncp)
  {

    zi=t(Xi)%*%yc #For cov(t(X)y,w)=0, w must be colinear to t(X)y ==> Z=t(X)y
    zi=as.vector(zi)

    Xsvd=svd(t(Xi)%*%Xi)
    XtXmoins1=Xsvd$v%*%diag(1/Xsvd$d)%*%t(Xsvd$u)
    #XtXmoins1=solve(t(Xi)%*%Xi)

    wLS=XtXmoins1%*%zi


    #Optimizing nu
    wLSs=sort(abs(wLS))

    wsp=(1:p)/p
    iz=which.min(abs(wsp-ppnu))

    delta=sign(XtXmoins1%*%zi)


    ###########
    nu=wLSs[iz]
    ###########

    # finding lambda, mu, given nu
    znu=wLS
    znu=sapply(znu,function(u) sign(u)*max(abs(u)-nu,0))

    # calculating w,t at the optimum
    w=znu

    #Finding t
    t=Xi%*%w
    t=t/d.spls.norm2(t)

    WW[,ic]=w
    TT[,ic]=t

    #Deflation
    Xi=Xi-t%*%t(t)%*%Xi

    #Coefficient vectors
    R=t(TT[,1:ic,drop=FALSE])%*%Xc%*%WW[,1:ic,drop=FALSE]
    R[row(R)>col(R)]<-0 # inserted for numerical stability

    L=backsolve(R,diag(ic))
    Bhat[,ic]=WW[,1:ic,drop=FALSE]%*%(L%*%(t(TT[,1:ic,drop=FALSE])%*%yc))

    intercept[ic] = ym - Xm %*% Bhat[,ic]
    zerovar[ic]=sum(Bhat[,ic]==0)

    #Predictions
    YY[,ic]=X %*% Bhat[,ic] + intercept[ic]
    RES[,ic]=y-YY[,ic]

    # results iteration
    if (verbose){
      cat('Dual PLS LS, ic=',ic,
          'nu=',nu,
          'nbzeros=',zerovar[ic], '\n')
    }
  }

  return(list(Xmean=Xm,scores=TT,loadings=WW,Bhat=Bhat,intercept=intercept,
              fitted.values=YY,residuals=RES,
              zerovar=zerovar))
}
