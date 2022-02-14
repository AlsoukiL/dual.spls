#' Dual Sparse Partial Least Squares (Dual-SPLS) regression for the least squares norm
#' @description
#' The function \code{d.spls.LS} performs dimensional reduction as in PLS1 methodology combined to variable selection via the
#' Dual-SPLS algorithm with the norm \deqn{\Omega(w)=\lambda \|N_1w\|_1 + \|Xw\|_2.}{\Omega(w)=\lambda ||N1w||_1 + ||Xw||_2.}
#' @usage d.spls.LS(X,y,ncp,ppnu,verbose=TRUE)
#' @param X a numeric matrix of predictors values of dimension \code{(n,p)}. Each row represents one observation and each column one predictor variable.
#' @param y a numeric vector or a one column matrix of responses. It represents the response variable for each observation.
#' @param ncp a positive integer. \code{ncp} is the number of Dual-SPLS components.
#' @param ppnu a positive real value, in \eqn{[0,1]}. \code{ppnu} is the desired
#' proportion of variables to shrink to zero for each component (see Dual-SPLS methodology).
#' @param verbose a Boolean value indicating whether or not to display the iterations steps. Default value is \code{TRUE}.
#' @details
#' The resulting solution for \eqn{w} and hence for the coefficients vector, in the case of \code{d.spls.LS}, has
#' a simple closed form expression (ref) deriving from the fact that \eqn{w} is collinear to a vector \eqn{z_{\nu}}{z.\nu} of coordinates
#' \deqn{z_{\nu_j}=\textrm{sign}(\hat{\beta}_{LS_j})(|\hat{\beta}_{LS_j}|-\nu)_+.}{z.\nu_j=sign(hat\beta.LS_j)(|hat\beta.LS_j|-\nu)_+.}
#' Here \eqn{\nu} is the threshold for which \code{ppnu} of
#' the absolute values of the coordinates of \eqn{\hat{\beta}_{LS}}{hat\beta.LS} are greater than \eqn{\nu} where \eqn{\hat{\beta}_{LS}=(X^TX)^{-1}X^Ty}{hat\beta.LS=inv(X^TX)X^Ty}.
#' Therefore, the least squares norm is only adapted to the situation where XT X is invertible.
#' At each step the singularity of XT X is tested by computing its condition number.
#' A finite large ratio means that the matrix is close to being singular.
#'
#' \eqn{N_1}{N1} does not intervene in the resolution step. Whoever, \eqn{z_{\nu_j}}{z.\nu_j} is constructed according \eqn{N_1}{N1}.
#' Therefore, proving that \eqn{N_1}{N1} exists is enough. (see references for more details.)
#' If the singularity condition is not verified, one might choose to apply the Dual-SPLS for the ridge norm.
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
#' \item{zerovar}{the vector of length \code{ncp} representing the number of variables shrank to zero per component.}
#' \item{ind_diff0}{the list of \code{ncp} elements representing the index of the none null regression coefficients elements.}
#' \item{type}{a character specifying the Dual-SPLS norm used. In this case it is \code{LS}. }
#' @author Louna Alsouki François Wahl
#' @seealso  [dual.spls::d.spls.ridge]
#'
#' @examples
#' ### load dual.spls library
#' library(dual.spls)
#' ### parameters
#' set.seed(14)
#' n <- 1000
#' p <- 40
#' nondes <- 100
#' sigmaondes <- 0.01
#' data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
#'
#' X <- data$X
#' y <- data$y
#'
#' #fitting the model
#' mod.dspls <- d.spls.LS(X=X,y=y,ncp=7,ppnu=0.9,verbose=TRUE)
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


d.spls.LS<- function(X,y,ncp,ppnu,verbose=TRUE)
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
  WW=matrix(0,p,ncp) # initialising W, the matrix of loadings
  TT=matrix(0,n,ncp) # initialising T, the matrix of scores
  Bhat=matrix(0,p,ncp) # initialising the matrix of coefficients
  YY=matrix(0,n,ncp) # initialising the matrix of estimated responses
  RES=matrix(0,n,ncp) # initialising the matrix of residues
  intercept=rep(0,ncp) # initialising intercept, the vector of intercepts
  zerovar=rep(0,ncp) # initialising zerovar, the vector of final number of zeros coefficients for each component
  listelambda=rep(0,ncp) # initialising listelambda, the vector of values of lambda
  ind.diff0=vector(mode = "list", length = ncp) # initializing ind0, the list of the index of the none zero coefficients
  names(ind.diff0)=paste0("in.diff0_", 1:ncp)

  ###################################
  # Dual-SPLS
  ###################################

  # each step ic in -for loop- determine the icth column or element of each element initialized
  Xi=Xc # initialising X for Deflation Step
  for (ic in 1:ncp)
  {

    zi=t(Xi)%*%yc
    zi=as.vector(zi)

    # computing the inverse of t(X)%*%X
    Xsvd=svd(Xi)
    len=min(p,n)
    #Xsvd=svd(t(Xi)%*%Xi)
    invD=1/Xsvd$d
    if (ic==1 && min(Xsvd$d)/max(Xsvd$d)<1e-16)
    {
      warning('XtX is close to being singular')
    }
    else
      if (ic>1 && min(Xsvd$d)/max(Xsvd$d[-((len-ic+2):len)])<1e-16)
      {
        warning('deflated XtX is close to being singular on component number ',ic )
        invD[((len-ic+2):len)]=0
      }
    XtXmoins1=Xsvd$v%*%diag(invD^2)%*%t(Xsvd$v)

    #Optimizing nu
    wLS=XtXmoins1%*%zi
    wLSs=sort(abs(wLS))
    wsp=(1:p)/p
    iz=which.min(abs(wsp-ppnu))
    delta=sign(XtXmoins1%*%zi)

    ###########
    nu=wLSs[iz]#
    ###########

    # finding lambda, mu, given nu
    znu=wLS
    znu=sapply(znu,function(u) sign(u)*max(abs(u)-nu,0))
    XZnu=Xi%*%znu
    XZnu2=d.spls.norm2(XZnu)
    #########
    mu=XZnu2 #
    ##############
    lambda=nu/mu #
    ##############

    # calculating w,t at the optimum
    w=znu
    WW[,ic]=w

    #Finding t
    t=Xi%*%w
    t=t/d.spls.norm2(t)
    TT[,ic]=t

    # deflation
    Xi=Xi-t%*%t(t)%*%Xi

    # coefficient vectors
    R=t(TT[,1:ic,drop=FALSE])%*%Xc%*%WW[,1:ic,drop=FALSE]
    R[row(R)>col(R)]<-0 # inserted for numerical stability

    L=backsolve(R,diag(ic))
    Bhat[,ic]=WW[,1:ic,drop=FALSE]%*%(L%*%(t(TT[,1:ic,drop=FALSE])%*%yc))

    # lambda
    listelambda[ic]=lambda

    # intercept
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
      cat('Dual PLS LS, ic=',ic,
          'nu=',nu,
          'nbzeros=',zerovar[ic], '\n')
    }
  }

  return(list(Xmean=Xm,scores=TT,loadings=WW,Bhat=Bhat,intercept=intercept,
              fitted.values=YY,residuals=RES,
              lambda=listelambda,zerovar=zerovar,ind.diff0=ind.diff0,type="LS"))
}
