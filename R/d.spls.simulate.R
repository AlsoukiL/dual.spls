#' Simulation of a data
#' @description
#' The function \code{d.spls.simulate} simulates \code{G} mixtures of \code{nondes} Gaussians from which it builds
#' a data set of predictors \code{X} and response \code{y} in a way that \code{X} can be divided into \code{G} groups and
#' the values of \code{y} depend on the values of \code{X}.
#' @usage d.spls.simulate(n=200,p=100,nondes=50,sigmaondes=0.05,sigmay=0.5)
#' @param n a positive integer. \code{n} is the number of observations. Default value is \code{200}.
#' @param p a numeric vector of length \code{G} representing the number of variables. Default value is \code{100}.
#' @param nondes a numeric vector of length \code{G}. \code{nondes} is the number of Guassians in each mixture. Default value is \code{50}.
#' @param sigmaondes a numeric vector of length \code{G}. \code{sigmaondes} is the standard deviation of the
#' Gaussians for each group \eqn{g \in [1,G]}. Default value is \code{0.05}.
#' @param sigmay a real value. \code{sigmay} is the uncertainty on \code{y}. Default value is \code{0.5}.
#' @details
#' The predictors matrix \code{X} is a concatunations of \code{G} predictors sub matrices. Each is computed using
#' a mixture of Gaussian i.e. summing the following Gaussians:
#' \deqn{A \exp{(-\frac{(\textrm{xech}-\mu)^2}{2 \sigma^2})}}.
#' Where
#' \itemize{
#' \item \eqn{A} is a numeric vector of random values between 0 and 1,
#' \item xech is a sequence of \eqn{p(g)} equally spaced values from 0 to 1. \eqn{p(g)} is the number
#' of variables of the sub matrix \eqn{g}, for \eqn{g \in \{1, \dots, G\}},
#' \item \eqn{\mu} is a random value in \eqn{[0,1]} representing the mean of the Gaussians,
#' \item \eqn{\sigma} is a positive real value specified by the user and representing the standard
#' deviation of the Gaussians.
#' }
#' The response vector \code{y} is a linear combination of the predictors to which we add a noise of uncertainty \code{sigmay}. It is computed as follows:
#'
#' \deqn{y_i=\textrm{sigmay} \times V_i +\sum\limits_{g=1}^G \sum\limits_{j=1}^4 j \times X^g_{ij}}
#' Where
#' \itemize{
#' \item \eqn{G} is the number of predictor sub matrices,
#' \item \eqn{i} is the index of the observation,
#' \item \eqn{X^g_{ij}} is the matrix of the columns \eqn{\textrm{init}_j, \dots, \textrm{fin}_j} from the \eqn{g} sub matrix \eqn{X^g},
#' \item \eqn{V} is a normally distributed vector of 0 mean and unitary standard deviation.
#' }
#' Noting that
#' \itemize{
#' \item \eqn{\textrm{init}_1= 0.2 \times p(g)} and \eqn{\textrm{fin}_1= 0.3 \times p(g)}
#' \item \eqn{\textrm{init}_2= \textrm{init}_1 \times 0.4 \times p(g)} and \eqn{\textrm{fin}_2= \textrm{fin}_1 \times 0.5 \times p(g)}
#' \item \eqn{\textrm{init}_3= \textrm{init}_2 \times 0.6 \times p(g)} and \eqn{\textrm{fin}_3= \textrm{fin}_2 \times 0.7 \times p(g)}
#' \item \eqn{\textrm{init}_4= \textrm{init}_3 \times 0.8 \times p(g)} and \eqn{\textrm{fin}_4= \textrm{fin}_3 \times 0.9 \times p(g)}
#' }
#' @return A \code{list} of the following attributes
#' \item{X}{the concatunated predictors matrix.}
#' \item{y}{the response vector.}
#' \item{y0}{the reponse vector without noise \code{sigmay}.}
#' \item{sigmay}{the uncertainty on \code{y}.}
#' \item{sigmaondes}{the standard deviation of the Gaussians.}
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
#' X2 <- Xb[,(p[1]+1):(p[1]+p[2])]
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
  if (length(nondes) != length(p))
  {
    if (length(nondes)==1)
    {
      warning("nondes is only specified for the first mixture. Same value is considered for the rest")
      nondes=rep(nondes,length(p))
    }
    else
    {
      stop("dimensions of nondes and p differ")
    }
  }

  if (length(sigmaondes) != length(p))
  {
    if (length(sigmaondes)==1)
    {
      warning("sigmaondes is only specified for the first mixture. Same value is considered for the rest")
      sigmaondes=rep(sigmaondes,length(p))
    }
    else
    {
      stop("dimensions of sigmaondes and p differ")
    }
  }
  X=matrix(0,n,sum(p)) #initializing the predictors matrix X

  ampl=matrix(runif(nondes[1]*n),nrow=n, ncol=nondes[1])
  modes=runif(nondes[1])
  xech=seq(0,1,length.out=p[1]) #setting the p variables

  for (j in 1:n){

    for (io in 1:nondes[1]){
      X[j,1:p[1]]=X[j,1:p[1]]+ampl[j,io]*exp(-(xech-modes[io])^2/(2*sigmaondes[1]^2))
    }
  }

  if (length(p)>1)
  {
    for (i in 2:length(p))
    {
      ampl=matrix(runif(nondes[i]*n),nrow=n, ncol=nondes[i])
      modes=runif(nondes[i])
      xech=seq(0,1,length.out=p[i]) #setting the p variables

      for (j in 1:n){
        for (io in 1:nondes[i]){
          X[j,(sum(p[1:(i-1)])+1):sum(p[1:i])]=X[j,(sum(p[1:(i-1)])+1):sum(p[1:i])]+ampl[j,io]*exp(-(xech-modes[io])^2/(2*sigmaondes[i]^2))
        }
      }
    }
  }

  y0=rep(0,n) #initializing the response vector without noise y0
  #Setting the interval limits for y0
  p1i=round(20/100*p[1])
  p1f=round(30/100*p[1])
  p2i=round(40/100*p[1])
  p2f=round(50/100*p[1])
  p3i=round(60/100*p[1])
  p3f=round(70/100*p[1])
  p4i=round(80/100*p[1])
  p4f=round(90/100*p[1])

  #Computing y0 as a sum of intervals of X
  for (j in 1:n){

    y0[j]=sum(X[j,p1i:p1f])+
      2* sum(X[j,p2i:p2f])+
      3* sum(X[j,p3i:p3f])+
      4* sum(X[j,p4i:p4f])
  }

  if (length(p)>1)
  {
    for (i in 2:length(p))
    {
      p1i=round(20/100*p[i])
      p1f=round(30/100*p[i])
      p2i=round(40/100*p[i])
      p2f=round(50/100*p[i])
      p3i=round(60/100*p[i])
      p3f=round(70/100*p[i])
      p4i=round(80/100*p[i])
      p4f=round(90/100*p[i])
      for (j in 1:n){

        y0[j]=y0[j]+sum(X[j,(sum(p[1:(i-1)])+p1i):(sum(p[1:(i-1)])+p1f)])+
          2* sum(X[j,(sum(p[1:(i-1)])+p2i):(sum(p[1:(i-1)])+p2f)])+
          3* sum(X[j,(sum(p[1:(i-1)])+p3i):(sum(p[1:(i-1)])+p3f)])+
          4* sum(X[j,(sum(p[1:(i-1)])+p4i):(sum(p[1:(i-1)])+p4f)])
      }
    }
  }
  #adding noise to y0
  y=y0+sigmay*rnorm(n)
  G=length(p)

  return(list(X=X,y=y,y0=y0,sigmay=sigmay,sigmaondes=sigmaondes,G=G))
}
