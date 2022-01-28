#' Dual Sparse Partial Least Squares (Dual-SPLS) regression for the lasso norm
#' @description
#' The function \code{d.spls.lasso} performs dimensional reduction as in the PLS1 methodology combined with variable selection via the
#' Dual-SPLS algorithm with the norm \deqn{\Omega(w)=\lambda \|w\|_1 + \|w\|_2.}{\Omega(w)=\lambda ||w||_1 + ||w||_2.}
#' @usage d.spls.lasso(X,y,ncp,ppnu,verbose=TRUE)
#' @param X a numeric matrix of predictors values of dimension \code{(n,p)}. Each row represents one observation and each column one predictor variable.
#' @param y a numeric vector or a one column matrix of responses. It represents the response variable for each observation.
#' @param ncp a positive integer. \code{ncp} is the number of Dual-SPLS components.
#' @param ppnu a positive real value, in \eqn{[0,1]}. \code{ppnu} is the desired
#' proportion of variables to shrink to zero for each component (see Dual-SPLS methodology).
#' @param verbose a Boolean value indicating whether or not to display the iterations steps. Default value is \code{TRUE}.
#' @details
#' The resulting solution for \eqn{w} and hence for the coefficients vector, in the case of \code{d.spls.lasso}, has
#' a simple closed form expression (ref) deriving from the fact that \eqn{w} is collinear to a vector \eqn{z_{\nu}}{z.{\nu}} of coordinates
#' \deqn{z_{\nu_j}=\textrm{sign}({z_j})(|z_j|-\nu)_+.}{z.\nu_j=sign(z_j)(|z_j|-\nu)_+.}
#' Here \eqn{\nu} is the threshold for which \code{ppnu} of
#' the absolute values of the coordinates of \eqn{z=X^Ty} are greater than \eqn{\nu}.
#'
#' @return A \code{list} of the following attributes
#' \item{Xmean}{the mean vector of the predictors matrix \code{X}.}
#' \item{scores}{the matrix of dimension \code{(n,ncp)} where \code{n} is the number of observations. The \code{scores} represents
#' the observations in the new component basis computed by the compression step
#' of the Dual-SPLS.}
#' \item{loadings}{the matrix of dimension \code{(p,ncp)} that represents the Dual-SPLS components.}
#' \item{Bhat}{the matrix of dimension \code{(p,ncp)} that regroups the regression coefficients for each component.}
#' \item{intercept}{the vector of intercept values for each component.}
#' \item{fitted.values}{the matrix of dimension \code{(n,ncp)} that represents the predicted values of \code{y}}
#' \item{residuals}{the matrix of dimension \code{(n,ncp)} that represents the residuals corresponding
#' to the difference between the responses and the fitted values.}
#' \item{lambda}{the vector of length \code{ncp} collecting the parameters of sparsity  used to fit the model at each iteration.}
#' \item{zerovar}{the vector of length \code{ncp} representing the number of variables shrank to zero per component.}
#' \item{ind_diff0}{the list of \code{ncp} elements representing the index of the none null regression coefficients elements.}
#' \item{type}{a character specifying the Dual-SPLS norm used. In this case it is \code{lasso}. }
#' @author Louna Alsouki François Wahl
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
#' main="Observed VS Predicted for 6 components")
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


d.spls.lasso<- function(X,y,ncp,ppnu,verbose=TRUE)
{

  ###################################
  # Dimensions
  ###################################
  n=length(y) # number of observations
  p=dim(X)[2] # number of variables

  ###################################
  # Centering Data
  ###################################
  Xm = apply(X, 2, mean) # mean of X
  Xc=X - rep(1,n) %*% t(Xm) # centering predictor matrix

  ym=mean(y) # mean of y
  yc=y-ym # centering response vector

  ###################################
  # Initialization
  ###################################
  WW=matrix(0,p,ncp) # initialising WW, the matrix of loadings
  TT=matrix(0,n,ncp) # initialising TT, the matrix of scores
  Bhat=matrix(0,p,ncp) # initialising Bhat, the matrix of coefficients
  YY=matrix(0,n,ncp) # initialising YY, the matrix of coefficients
  RES=matrix(0,n,ncp) # initialising RES, the matrix of coefficients
  intercept=rep(0,ncp) # initialising intercept, the vector of intercepts
  zerovar=rep(0,ncp) # initialising zerovar, the vector of final number of zeros coefficients for each component
  listelambda=rep(0,ncp) # initialising listelambda, the vector of values of lambda
  ind.diff0=vector(mode = "list", length = ncp) # initializing ind0, the list of the index of the none zero coefficients
  names(ind.diff0)=paste0("in.diff0_", 1:ncp)

  ###################################
  # Dual-SPLS
  ###################################

  # each step ic in -for loop- determine the icth column or element of each element initialized
  Xdef=Xc # initializing X for Deflation Step
  for (ic in 1:ncp)
  {

    Z=t(Xdef)%*%yc
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

    # finding W
    w=(mu/(nu*Znu1 + mu^2))*Znu
    WW[,ic]=w

    # finding T
    t=Xdef%*%w
    t=t/d.spls.norm2(t)
    TT[,ic]=t

    # deflation
    Xdef=Xdef-t%*%t(t)%*%Xdef

    # coefficient vectors
    R=t(TT[,1:ic,drop=FALSE])%*%Xc%*%WW[,1:ic,drop=FALSE]
    R[row(R)>col(R)]<-0 # inserted for numerical stability

    L=backsolve(R,diag(ic))
    Bhat[,ic]=WW[,1:ic,drop=FALSE]%*%(L%*%(t(TT[,1:ic,drop=FALSE])%*%yc))

    # lambda
    listelambda[ic]=lambda

    # intercept
    intercept[ic] = ym - Xm %*% Bhat[,ic]

    # zerovar
    zerovar[ic]=sum(Bhat[,ic]==0)

    # ind.diff0
    name=paste("in.diff0_",ic,sep="")
    ind.diff0[name]=list(which(Bhat[,ic]!=0))

    # predictions
    YY[,ic]=X %*% Bhat[,ic] + intercept[ic]

    # residuals
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
              lambda=listelambda,zerovar=zerovar,ind.diff0=ind.diff0,type="lasso"))
}
