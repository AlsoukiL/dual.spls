#' Plots the coefficient curve of a Dual-SPLS regression
#' @description
#' The function \code{dual.spls.plot} provides the regression coefficient curves of a Dual-SPLS model for a specified number of components
#' and the mean of the original data plot.
#' @usage d.spls.plot(mod.dspls,ncomp)
#' @param mod.dspls is a fitted Dual-SPLS object.
#' @param ncomp a positive integer or a numeric vector of the number of Dual-SPLS components to consider.
#' @details
#' The plots allow the visualization of the results by comparing the mean of the original data to the coefficients regression in the
#' case of a Dual-SPLS regression. The plots provided correspond to the Dual-SPLS coefficients for each \code{ncomp} desired.
#'
#' @return no return value
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
#' #fitting the Dual-SPLS lasso model
#'
#' mod.dspls.lasso <- d.spls.lasso(X=X,y=y,ncp=10,ppnu=0.9,verbose=TRUE)
#'
#' ncomp=c(5,6,7)
#' d.spls.plot(mod.dspls.lasso,ncomp)
#' @export

d.spls.plot<-function (mod.dspls,ncomp)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  # mean of the original data
  Xmean=mod.dspls$Xmean
  p=length(Xmean)

  par(mfrow=c(2,1))
  # raw data plot
  plot(Xmean,type='l',main=paste("Mean of the original data"), ylab='X',xlab='index',col=1)

  # Dual-SPLS plot
  for (i in 1:length(ncomp))
  {
  nz=mod.dspls$zerovar[ncomp[i]]
  plot(1:p,mod.dspls$Bhat[,ncomp[i]],type='l',
       main=paste(" Dual-SPLS (",mod.dspls$type,") for ncp = ", ncomp[i], " #0coef =", nz, "/", p),
       ylab='coefficients',xlab='index' )
  inonz=which(mod.dspls$Bhat[,ncomp[i]]!=0)
  points(inonz,mod.dspls$Bhat[inonz,ncomp[i]],col='red',pch=19,cex=0.5)
  legend("topright", legend ="non null values", bty = "n", cex = 0.8, col = "red",pch=19)
  }
}
