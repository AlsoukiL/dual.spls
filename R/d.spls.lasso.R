#' Dual Sparse Partial Least Squares (Dual-SPLS) regression for the norm \eqn{\Omega(w)=\lambda \|w\|_1 + \|w\|_2}
#' @description
#' The function \code{d.spls.lasso} performs dimensional reduction as in PLS methodology combined to variable selection via the
#' Dual-SPLS algorithm with the norm \eqn{\Omega(w)=\lambda \|w\|_1 + \|w\|_2}.
#' @usage d.spls.lasso(X,y,ncp,ppnu,verbose=FALSE)
#' @param X a numeric matrix of predictors values. Each row represents one observation and each column a predictor variable.
#' @param y a numeric vector or a one column matrix of responses. It represents the response variable for each observation.
#' @param ncp a positive integer. \code{ncp} is the number of Dual-SPLS components.
#' @param ppnu a positive real value, in \code{[0,1]}. \code{ppnu} is the desired
#' proportion of variables to shrink to zero for each component (see Dual-SPLS methodology).
#' @param verbose a boolean value indicating whether or not to diplay the iterations steps.
#' @details
#' This procedure computes latent sparse components that are used in a regression model.
#' The optimization problem of the Dual-SPLS regression for the norm \eqn{\Omega(w)=\lambda \|w\|_1 + \|w\|_2}
#' comes from the Dual norm defintion of \eqn{\Omega(w)}. Indeed, we are searching for \code{w} that goes
#' with
#' \deqn{\Omega^*(z)=max_w(z^Tw) \text{ s.t. } \Omega(w)=1}
#' Noting that \eqn{\lambda} is the initial shrinkage parameter that imposes sparsity, the Dual-SPLS does not rely
#' on the value of \eqn{\lambda} but instead proceeds adaptively by choosing the propotion of zeros that the user
#' would like to impose in the coefficients. Which leads to the compuation of a secondary shrinkage parameter \eqn{\nu}.
#'
#' The solution of the optimization problem is
#' \deqn{\dfrac{w_j}{\|w\|_2}=\dfrac{1}{\mu} \delta_j (|z_j|-\nu)_+}
#' Where
#' \itemize{
#' \item \eqn{\delta_j} is the vector of signs of \code{z} and \code{w}
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
#' \item{lambda}{the vector of length \eqn{ncp} collecting the parameters of sparsity  used to fit the model at each iteration.}
#' \item{zerovar}{the vector of length \eqn{ncp} representing the number of variables shrinked to zero per component.}
#' @author Louna Alsouki François Wahl
#' @seealso `browseVignettes("dual.spls")`
#'
#' @examples
#' ### load dual.spls library
#' library(dual.spls)
#' ### constructing the simulated example
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
#' mod.dspls <- d.spls.lasso(X=X,y=y,ncp=10,ppnu=0.9,verbose=TRUE)
#'
#' str(mod.dspls)
#'
#' ### plotting the observed values VS predicted values for 6 components
#' plot(y,mod.dspls$fitted.values[,6], xlab="Observed values", ylab="Predicted values",
#'  main="Observed VS Predicted for 6 components")
#' points(-1000:1000,-1000:1000,type='l')
#'
#' ### plotting the regression coefficients
#' par(mfrow=c(3,1))
#'
#' i=6
#' nz=mod.dspls$zerovar[i]
#' plot(1:dim(X)[2],mod.dspls$Bhat[,i],type='l',
#'     main=paste(" Dual-SPLS (lasso), ncp =", i, " #0coef =", nz, "/", dim(X)[2]),
#'     ylab='',xlab='' )
#' inonz=which(mod.dspls$Bhat[,i]!=0)
#' points(inonz,mod.dspls$Bhat[inonz,i],col='red',pch=19,cex=0.5)
#' legend("topright", legend ="non null values", bty = "n", cex = 0.8, col = "red",pch=19)
#' @export


d.spls.lasso<- function(X,y,ncp,ppnu,verbose=FALSE)
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
  listelambda=rep(0,ncp)
  ###################################
  # Dual-SPLS
  ###################################

  #Each step ic in -for loop- determine the icth column of each W, T and Bhat
  Xdef=Xc #Initialising X for Deflation Step
  for (ic in 1:ncp)
  {

    Z=t(Xdef)%*%yc #For cov(t(X)y,w)=0, w must be colinear to t(X)y ==> Z=t(X)y
    Z=as.vector(Z)

    #Optimizing nu
    Zs=sort(abs(Z))
    Zsp=(1:p)/p
    iz=which.min(abs(Zsp-ppnu))
    ###########
    nu=Zs[iz] #
    ###########

    # finding lambda, mu, given nu
    Znu=sapply(Z,function(u) sign(u)*max(abs(u)-nu,0))
    Znu2=d.spls.norm2(Znu)
    Znu1=d.spls.norm1(Znu)
    #########
    mu=Znu2 #
    ##############
    lambda=nu/mu #
    ##############

    # calculating w,t at the optimum
    w=(mu/(nu*Znu1 + mu^2))*Znu

    #Finding T
    t=Xdef%*%w
    t=t/d.spls.norm2(t)

    WW[,ic]=w
    TT[,ic]=t

    #Deflation
    Xdef=Xdef-t%*%t(t)%*%Xdef

    #Coefficient vectors
    R=t(TT[,1:ic,drop=FALSE])%*%Xc%*%WW[,1:ic,drop=FALSE]
    R[row(R)>col(R)]<-0 # inserted for numerical stability

    L=backsolve(R,diag(ic))
    Bhat[,ic]=WW[,1:ic,drop=FALSE]%*%(L%*%(t(TT[,1:ic,drop=FALSE])%*%yc))

    listelambda[ic]=lambda
    intercept[ic] = ym - Xm %*% Bhat[,ic]
    zerovar[ic]=sum(Bhat[,ic]==0)

    #Predictions
    YY[,ic]=X %*% Bhat[,ic] + intercept[ic]
    RES[,ic]=y-YY[,ic]

    # results iteration
    if (verbose){
      cat('Dual PLS ic=',ic,'lambda=',lambda,
          'mu=',mu,'nu=',nu,
          'nbzeros=',zerovar[ic], '\n')
    }
  }

  return(list(Xmean=Xm,scores=TT,loadings=WW,Bhat=Bhat,intercept=intercept,
              fitted.values=YY,residuals=RES,
              lambda=listelambda,zerovar=zerovar))
}
