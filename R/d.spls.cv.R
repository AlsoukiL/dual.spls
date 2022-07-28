#' Determination of the number of latent components to be used in a Dual-SPLS regression
#' @description
#' The function \code{d.spls.cv} uses the cross validation approach described in Boulesteix and Strimmer (2005) (see in references) in order to
#' choose the most adequat number of latent components for a Dual-SPLS regression.
#' @usage d.spls.cv(X,Y,ncomp,dspls="lasso",ppnu,nu2,nrepcv=30,pctcv=70,indG,gamma)
#' @param X a numeric matrix of predictors values of dimension \code{(n,p)}. Each row represents one observation and each column one predictor variable.
#' @param Y a numeric vector or a one column matrix of responses. It represents the response variable for each observation.
#' @param ncomp a positive integer or a numeric vector of the number of Dual-SPLS components to choose from.
#' @param dspls the norm type of the Dual-SPLS regression applied. Default value is \code{lasso}. Options are \code{pls}, \code{LS},
#' \code{ridge}, \code{GLA}, \code{GLB} and \code{GLC}.
#' @param ppnu a positive real value, in \eqn{[0,1]}. \code{ppnu} is the desired
#' proportion of variables to shrink to zero for each component (see Dual-SPLS methodology).
#' @param nu2 a positive real value. \code{nu2} is a constraint parameter used in the ridge norm.
#' @param nrepcv a positive integer indicating the number of cross-validation iterations to be performed. Default value is 30.
#' @param pctcv a positive real value in \eqn{[0,100]} representing the percentage of observation to be considered in for the
#' calibration set at each CV iteration. Default value is 70.
#' @param indG a numeric vector of group index for each observation. It is used in the cases of the group lasso norms.
#' @param gamma a numeric vector of the norm \eqn{\Omega} of each \eqn{w_g} in case \code{GLB}.
#' @details
#' The procedure is described in the Boulesteix and Strimmer. It is based on randomly selecting, \code{pctcv%} of calibration observations at each
#' cross validation iteration and performing the Dual-SPLS regression. The rest of the observation are used as a validation and the
#' errors are computed accordingly for each components. \code{nrepcv} iterations are done and the mean squared of each of the \code{nrepcv} errors for each
#' component are computed. The latent component with the smallest mean value is selected as the best.
#' @return A \code{integer} representing the best number of latent components to be used in the Dual-SPLS regression based on the cross validation procedure.
#' @references
#' A. L. Boulesteix and K. Strimmer (2005). Predicting Transcription Factor Activities from Combined Analysis of Microarray and ChIP Data: A Partial Least Squares Approach. \cr
#' \cr
#' H. Wold. Path Models with Latent Variables: The NIPALS Approach. In H.M. Blalock et al., editor, Quantitative Sociology: International Perspectives on Mathematical and Statistical Model Building, pages 307–357. Academic Press, 1975.
#' @author Louna Alsouki François Wahl
#'
#' @examples
#' ### load dual.spls library
#' library(dual.spls)
#' ### constructing the simulated example
#' oldpar <- par(no.readonly = TRUE)
#' n <- 100
#' p <- 50
#' nondes <- 20
#' sigmaondes <- 0.5
#' data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
#'
#' X <- data$X
#' y <- data$y
#'
#' #fitting the PLS model
#' ncomp_PLS <- d.spls.cv(X=X,Y=y,ncomp=10,dspls="pls",nrepcv=20,pctcv=75)
#' mod.dspls.pls <- d.spls.pls(X=X,y=y,ncp=ncomp_PLS,verbose=TRUE)
#'
#' str(mod.dspls.pls)
#'
#' ### plotting the observed values VS predicted values for ncomp components
#' plot(y,mod.dspls.pls$fitted.values[,ncomp_PLS], xlab="Observed values", ylab="Predicted values",
#'  main=paste("Observed VS Predicted for ", ncomp_PLS," components"))
#' points(-1000:1000,-1000:1000,type='l')
#'
#' ### plotting the regression coefficients
#' par(mfrow=c(3,1))
#'
#' i=ncomp_PLS
#' plot(1:dim(X)[2],mod.dspls.pls$Bhat[,i],type='l',
#'     main=paste(" Dual-SPLS (PLS), ncp =", i,
#'     ylab='',xlab='' ))
#'
#'
#' #fitting the Dual-SPLS lasso model
#'
#' ncomplasso <- d.spls.cv(X=X,Y=y,ncomp=10,dspls="lasso",ppnu=0.9,nrepcv=20,pctcv=75)
#' mod.dspls.lasso <- d.spls.lasso(X=X,y=y,ncp=ncomplasso,ppnu=0.9,verbose=TRUE)
#'
#' str(mod.dspls.lasso)
#'
#' ### plotting the observed values VS predicted values for ncomp components
#' plot(y,mod.dspls.lasso$fitted.values[,ncomplasso], xlab="Observed values", ylab="Predicted values",
#' main=paste("Observed VS Predicted for ", ncomplasso," components"))
#' points(-1000:1000,-1000:1000,type='l')
#'
#' ### plotting the regression coefficients
#' par(mfrow=c(3,1))
#'
#' i=ncomplasso
#' nz=mod.dspls.lasso$zerovar[i]
#' plot(1:dim(X)[2],mod.dspls.lasso$Bhat[,i],type='l',
#'     main=paste(" Dual-SPLS (lasso), ncp =", i, " #0coef =", nz, "/", dim(X)[2]),
#'     ylab='',xlab='' )
#' inonz=which(mod.dspls.lasso$Bhat[,i]!=0)
#' points(inonz,mod.dspls.lasso$Bhat[inonz,i],col='red',pch=19,cex=0.5)
#' legend("topright", legend ="non null values", bty = "n", cex = 0.8, col = "red",pch=19)
#' par(oldpar)
#' @export

d.spls.cv<-function (X, Y, ncomp,dspls="lasso",ppnu,nu2, nrepcv = 30, pctcv = 70,indG,gamma)
{
  n=dim(X)[1]
  ncal=floor(n * pctcv/100)
  cvcal=matrix(0, ncal, nrepcv)
  if (length(ncomp) == 1) {
    if (ncomp == 1)
      return(ncomp)
    else ncomp=1:ncomp
  }
  for (i in 1:nrepcv) {
    cvcal[, i]=sample(n, ncal)
  }
  cvcal=as.data.frame(cvcal)
  errorcv=sapply(cvcal, FUN = d.spls.errorcv, X,Y, ncomp,dspls,ppnu,nu2,indG,gamma)
  errorcv=errorcv[ncomp,]
  meanerror=apply(errorcv, MARGIN = 1, FUN = mean)
  ncomp=ncomp[which.min(meanerror)]
  return(ncomp)
}
