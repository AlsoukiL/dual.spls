#' Dual Sparse Partial Least Squares (Dual-SPLS) regression for the ridge norm
#' @description
#' The function \code{d.spls.lasso} performs dimensional reduction as in PLS methodology combined to variable selection via the
#' Dual-SPLS algorithm with the norm \deqn{\Omega(w)=\lambda_1 \|w\|_1 +\lambda_2 \|Xw\|_2 + \|w\|_2.}{\Omega(w)=\lambda1 ||w||_1 +\lambda2 ||Xw||_2 + ||w||_2.}.
#' In the algorithm, the parameters \eqn{\lambda}, \eqn{\lambda_1}{\lambda1} and \eqn{\lambda_2}{\lambda2} are transformed into more meaningful values, \code{ppnu} and \eqn{\nu_2}{\nu2}.
#' @usage d.spls.ridge(X,y,ncp,ppnu,nu2,verbose=TRUE)
#' @param X a numeric matrix of predictors values of dimension \code{(n,p)}. Each row represents one observation and each column one predictor variable.
#' @param y a numeric vector or a one column matrix of responses. It represents the response variable for each observation.
#' @param ncp a positive integer. \code{ncp} is the number of Dual-SPLS components.
#' @param ppnu a positive real value, in \eqn{[0,1]}. \code{ppnu} is the desired
#' proportion of variables to shrink to zero for each component (see Dual-SPLS methodology).
#' @param nu2 a positive real value. \code{nu2} is a regularization parameter on \eqn{X^TX}{XTX}.
#' @param verbose a Boolean value indicating whether or not to display the iterations steps. Default value is \code{TRUE}.
#' @details
#' The resulting solution for \eqn{w} and hence for the coefficients vector, in the case of \code{d.spls.ridge}, has
#' a simple closed form expression (ref) deriving from the fact that \eqn{w} is collinear to a vector \eqn{z_{\nu_1}}{z.\nu1} of coordinates
#' \deqn{z_{\nu_1,j}=\textrm{sign}(z_{X,\nu_2,j})(|z_{X,\nu_2,j}|-\nu_1)_+.}{z.\nu1_j=sign(z.X.\nu2_j)(|z.X.\nu2_j|-\nu1)_+.}
#' Here \eqn{\nu_1}{\nu1} is the threshold for which \code{ppnu} of
#' the absolute values of the coordinates of \eqn{z_{X,\nu_2}}{z.X.\nu2} are greater than \eqn{\nu_1}{\nu1} and \eqn{z_{X,\nu_2}=(\nu_2 X^TX + I_p)^{-1}X^Ty}{z.X.\nu2=inv(\nu2 X^TX + I) X^Ty}.
#' Therefore, the ridge norm is beneficial to the situation where \eqn{X^TX}{XTX} is singular. If \eqn{X^TX}{XTX} is invertible, one
#' can choose to use the Dual-SPLS for the least squares norm instead.
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
#'  to the difference between the responses and the fitted values.}
#' \item{lambda1}{the vector of length \code{ncp} collecting the parameters of sparsity used to fit the model at each iteration.}
#' \item{zerovar}{the vector of length \code{ncp} representing the number of variables shrank to zero per component.}
#' \item{ind_diff0}{the list of \code{ncp} elements representing the index of the none null regression coefficients elements.}
#' \item{type}{a character specifying the Dual-SPLS norm used. In this case it is \code{ridge}. }
#' @author Louna Alsouki François Wahl
#' @seealso [dual.spls::d.spls.LS]
#'
#' @examples
#' ### load dual.spls library
#' library(dual.spls)
#' ### parameters
#' oldpar <- par(no.readonly = TRUE)
#' n <- 200
#' p <- 100
#' nondes <- 150
#' sigmaondes <- 0.01
#' data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
#'
#' X <- data$X
#' y <- data$y
#'
#'
#' #fitting the model
#' mod.dspls <- d.spls.ridge(X=X,y=y,ncp=10,ppnu=0.9,nu2=100,verbose=TRUE)
#'
#' str(mod.dspls)
#'
#' ### plotting the observed values VS predicted values
#' plot(y,mod.dspls$fitted.values[,6], xlab="Observed values", ylab="Predicted values",
#' main="Observed VS Predicted for 6 components")
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
#' par(oldpar)
#' @export


d.spls.ridge<- function(X,y,ncp,ppnu,nu2,verbose=TRUE)
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
  # Initialisation
  ###################################
  WW=matrix(0,p,ncp) # initialising WW, the matrix of loadings
  TT=matrix(0,n,ncp) # initialising TT, the matrix of scores
  Bhat=matrix(0,p,ncp) # initialising Bhat, the matrix of coefficients
  YY=matrix(0,n,ncp) # initialising YY, the matrix of coefficients
  RES=matrix(0,n,ncp) # initialising RES, the matrix of coefficients
  intercept=rep(0,ncp) # initialising intercept, the vector of intercepts
  zerovar=rep(0,ncp) # initialising zerovar, the vector of final number of zeros coefficients for each component
  listelambda1=rep(0,ncp) # initialising listelambda1, the vector of values of lambda1
  ind.diff0=vector(mode = "list", length = ncp) # initializing ind0, the list of the index of the none zero coefficients
  names(ind.diff0)=paste0("in.diff0_", 1:ncp)

  ###################################
  # Dual-SPLS
  ###################################

  # each step ic in -for loop- determine the icth column or element of each element initialized
  Xdef=Xc # initialising X for Deflation Step
  for (ic in 1:ncp)
  {

    # Z=X^Ty
    z=t(Xdef)%*%yc
    z=as.vector(z)

    # computing z12
    # computing (nu2 N2T N2 + I)^(-1)
    N2=Xdef
    N2TN2=t(Xdef)%*%Xdef
    temp=nu2*N2TN2+diag(p)
    Xsvd=svd(temp)
    inv=Xsvd$v%*%diag(1/Xsvd$d)%*%t(Xsvd$u)

    # computing delta
    delta=sign(z)

    ZN=inv%*%z

    # optimizing nu
    zs=sort(abs(ZN))
    zsp=(1:p)/p
    iz=which.min(abs(zsp-ppnu))
    ###########
    nu1=zs[iz] #
    # ###########

    z12=sapply(ZN,function(u) sign(u)*max(abs(u)-nu1,0))

    # finding lambda1, mu, given nu
    # computing mu
    #########
    mu=d.spls.norm2(z12) #
    ##############
    lambda1=nu1/mu #
    ##############

    # calculating w,t at the optimum
    w=(mu/(nu2*d.spls.norm2(N2%*%z12)^2+nu1*d.spls.norm1(z12) + mu^2))*z12
    WW[,ic]=w

    #Finding T
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

    #lambda1
    listelambda1[ic]=lambda1

    #intercept
    intercept[ic] = ym - Xm %*% Bhat[,ic]

    # predictions
    YY[,ic]=X %*% Bhat[,ic] + intercept[ic]

    # residuals
    RES[,ic]=y-YY[,ic]

    #zerovar
    zerovar[ic]=sum(Bhat[,ic]==0)

    # ind.diff0
    name=paste("in.diff0_",ic,sep="")
    ind.diff0[name]=list(which(Bhat[,ic]!=0))

    # results iteration
    if (verbose){
      cat('Dual PLS ic=',ic,'lambda1=',lambda1,
          'mu=',mu,'nu2=',nu2,
          'nbzeros=',zerovar[ic], '\n')
    }

  }

  return(list(Xmean=Xm,scores=TT,loadings=WW,Bhat=Bhat,intercept=intercept,
              fitted.values=YY,residuals=RES,
              lambda1=listelambda1,zerovar=zerovar,ind.diff0=ind.diff0,type="ridge"))

}
