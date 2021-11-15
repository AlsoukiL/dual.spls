#' Simulation of a data
#' @description
#' The function \code{d.spls.simulate} simulates one or two mixtures of \code{nondes} Gaussians representing one or two
#' predictors matrices \code{X} or \code{X1} and \code{X2}
#' of \code{n} observations and \code{p} or \code{p1} and \code{p2} variables and a response vector \code{y}.
#' @usage d.spls.simulate(n=200,p=100,nondes=50,sigmaondes=0.05,sigmay=0.5)
#' @param n a positive integer. \code{n} is the number of observations. Default value is 200.
#' @param p a positive real value or a vector of length 2 representing the number of variables. Default value is 100.
#' @param nondes a positive integer or a vector of length 2. \code{nondes} is the number of Guassians in each mixture. Default value is 50.
#' @param sigmaondes a positive real value or a vector of length 2. \code{sigmaondes} is the standard deviation of the Gaussians. Default value is 0.05.
#' @param sigmay a real value. \code{sigmay} is the uncertainty on \code{y}. Default value is 0.5.
#' @details
#' The predictors matrices \code{X} or \code{X1} and \code{X2} represent each a mixture of \code{nondes} Gaussians characterized by means and
#' amplitudes chosen randomly and a standard deviation set to \code{sigmaondes}.
#'
#' The response vector \code{y} is the sum of 4 intervals of each predictor matrix to which we add noise of \code{sigmay}.
#' @return A \code{list} of the following attributes
#' \item{X}{the preditors matrix. If length(p)=2, it is the concatunated predictor matrix.}
#' \item{y}{the response vector.}
#' \item{y0}{the reponse vector without noise \code{sigmay}.}
#' \item{sigmay}{the uncertainty on \code{y}.}
#' \item{sigmaondes}{the standard deviation of the Gaussians.}
#' }
#' @author Louna Alsouki François Wahl
#' @seealso `browseVignettes("dual.spls")`
#'
#' @examples
#' ### load dual.spls library
#' library(dual.spls)
#' ####one predictors matrix
#' ### parameters
#' n <- 100
#' p <- 50
#' nondes <- 20
#' sigmaondes <- 0.5
#' data1=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
#'
#' Xa <- data1$X
#' ya <- data1$y
#'
#' ###plotting the data
#' plot(Xa[1,],type='l',ylim=c(0,max(Xa)),main='Data', ylab='Xa',col=1)
#' for (i in 2:n){
#'  lines(Xa[i,],col=i)
#'}
#'
#'####two predictors matrix
#' ### parameters
#' n <- 100
#' p <- c(50,100)
#' nondes <- c(20,30)
#' sigmaondes <- c(0.05,0.02)
#' data2=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
#'
#' Xb <- data2$X
#' X1 <- Xb[,(1:p[1])]
#' X2 <- Xb[,(p[1]+1):p[2]]
#' yb <- data2$y
#'
#' ###plotting the data
#' plot(Xb[1,],type='l',ylim=c(0,max(Xb)),main='Data', ylab='Xb',col=1)
#' for (i in 2:n){
#'  lines(Xb[i,],col=i)
#'}
#'
#' ###plotting the data
#' plot(X1[1,],type='l',ylim=c(0,max(X1)),main='Data X1', ylab='X1',col=1)
#' for (i in 2:n){
#'  lines(X1[i,],col=i)
#'}
#'
#' ###plotting the data
#' plot(X2[1,],type='l',ylim=c(0,max(X2)),main='Data X2', ylab='X2',col=1)
#' for (i in 2:n){
#'  lines(X2[i,],col=i)
#'}
#' @export


d.spls.simulate<- function(n=200,p=100,nondes=50,sigmaondes=0.05,sigmay=0.5)
{

  if (length(p)==1)
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

  }

  if (length(p)==2)
  {
    if (length(nondes)==1 || length(sigmaondes) ==1)
    {
      warning("nondes or sigmaondes is not specified for X2, same value is set for both mixtures")
      if (length(nondes) == 1) {nondes=c(nondes,nondes)}
      if (length(sigmaondes) == 1) {sigmaondes=c(sigmaondes,sigmaondes)}

    }
    p1=p[1]
    p2=p[2]
    nondes1=nondes[1]
    nondes2=nondes[2]
    sigmaondes1=sigmaondes[1]
    sigmaondes2=sigmaondes[2]

    ampl1=matrix(runif(nondes1*n),nrow=n, ncol=nondes1)
    modes1=runif(nondes1)

    ampl2=matrix(runif(nondes2*n),nrow=n, ncol=nondes2)
    modes2=runif(nondes2)

    X1=matrix(0,n,p1)
    X2=matrix(0,n,p2)
    y0=rep(0,n)
    xech1=seq(0,1,length.out=p1)
    xech2=seq(0,1,length.out=p2)


    for (i in 1:n){

      for (io in 1:nondes1){
        X1[i,]=X1[i,]+ampl1[i,io]*exp(-(xech1-modes1[io])^2/(2*sigmaondes1^2))
      }
    }

    for (i in 1:n){

      for (io in 1:nondes2){
        X2[i,]=X2[i,]+ampl2[i,io]*exp(-(xech2-modes2[io])^2/(2*sigmaondes2^2))
      }
    }

    p1i1=round(20/100*p1)
    p1f1=round(30/100*p1)
    p2i1=round(45/100*p1)
    p2f1=round(55/100*p1)
    p3i1=round(55/100*p1)
    p3f1=round(60/100*p1)
    p4i1=round(70/100*p1)
    p4f1=round(80/100*p1)

    p1i2=round(20/100*p2)
    p1f2=round(30/100*p2)
    p2i2=round(45/100*p2)
    p2f2=round(55/100*p2)
    p3i2=round(55/100*p2)
    p3f2=round(60/100*p2)
    p4i2=round(70/100*p2)
    p4f2=round(80/100*p2)

    for (i in 1:n){

      y0[i]=sum(X1[i,p1i1:p1f1])+
        2* sum(X1[i,p2i1:p2f1])+
        3* sum(X1[i,p3i1:p3f1])+
        4* sum(X1[i,p4i1:p4f1])

      y0[i]=y0[i]+sum(X2[i,p1i2:p1f2])+
        2* sum(X2[i,p2i2:p2f2])+
        3* sum(X2[i,p3i2:p3f2])+
        4* sum(X2[i,p4i2:p4f2])
    }

    X=cbind(X1,X2)
  }

  #adding noise to y0
  y=y0+sigmay*rnorm(n)

  return(list(X=X,y=y,y0=y0,sigmay=sigmay,sigmaondes=sigmaondes))
}
