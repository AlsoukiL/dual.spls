#' Simulation of a Benchmark data
#' @description
#' The function \code{BCHMK} simulates a mixture of \code{nondes} Gaussians representing a predictors matrix \code{X}
#' of \code{n} observations and \code{p} variables and a response vector \code{y}.
#' @usage
#' \code{BCHMK(n=200,p=100,nondes=50,sigmaondes=0.05,sigmay=0.5)}
#' @param n a positive integer. \code{n} is the number of observations.
#' @param p a positive integer. \code{p} is the number of variables.
#' @param nondes a positive integer. \code{p} is the number of Guassians in the mixture.
#' @param sigmaondes a positive real value. \code{sigmaondes} is the standard deviation of the Gaussians.
#' @param sigmay a real value. \code{sigmay} is the uncertainty on \code{y}.
#' @details
#' The predictors matrix \code{X} represents a mixture of \code{nondes} Gaussians characterized by means and
#' amplitudes chosen randomly and a standard deviation set to \code{sigmaondes}.
#'
#' The response vector \code{y} is the sum of 4 intervals of \code{X} to which we add noise of \code{sigmay}.
#' #' @return A \code{list} of the following attributes
#' @param X the preditors matrix.
#' @param y the response vector.
#' @param y0 the reponse vector without noise \code{sigmay}.
#' @param sigmay the uncertainty on \code{y}.
#' @param sigmaondes the standard deviation of the Gaussians.
#' @author François Wahl Louna Alsouki
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
#' data.benchmark=BCHMK(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
#'
#' X <- data.benchmark$X
#' y <- data.benchmark$y
#'
#' ###plotting the data
#' plot(X[1,],type='l',ylim=c(0,max(X)),main='Benchmark data', ylab='X',col=1)
#' for (i in 2:n){
#'  lines(X[i,],col=i)
#'}
#' @export


BCHMK<- function(n=200,p=100,nondes=50,sigmaondes=0.05,sigmay=0.5)
{

  ampl=matrix(runif(nondes*n),nrow=n, ncol=nondes) #computing random amplitudes
  modes=runif(nondes) #computing random means

  X=matrix(0,n,p) #initializing the predictors matrix X
  y0=rep(0,n) #initializing the response vector without noise y0
  xech=seq(0,1,length.out=p) #setting the p variables


  #combining the Guassian mixtures into X
  for (i in 1:n){

    for (io in 1:nondes){
      X[i,]=X[i,]+ampl[i,io]*exp(-(xech-modes[io])^2/(2*sigmaondes^2))
    }
  }

  #Setting the interval limits for y0
  p1i=round(20/100*p)
  p1f=round(30/100*p)
  p2i=round(45/100*p)
  p2f=round(55/100*p)
  p3i=round(55/100*p)
  p3f=round(60/100*p)
  p4i=round(70/100*p)
  p4f=round(80/100*p)

  #Computing y0 as a sum of intervals of X
  for (i in 1:n){

    y0[i]=sum(X[i,p1i:p1f])+
      2* sum(X[i,p2i:p2f])+
      3* sum(X[i,p3i:p3f])+
      4* sum(X[i,p4i:p4f])
  }

  #adding noise to y0
  y=y0+sigmay*rnorm(n)

  return(list(X=X,y=y,y0=y0,sigmay=sigmay,sigmaondes=sigmaondes))
}
