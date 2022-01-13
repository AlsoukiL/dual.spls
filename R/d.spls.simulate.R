#' Simulation of a data
#' @description
#' The function \code{d.spls.simulate} simulates \code{G} mixtures of \code{nondes} Gaussians from which it builds
#' a data set of predictors \code{X} and response \code{y} in a way that \code{X} can be divided into \code{G} groups and
#' the values of \code{y} depend on the values of \code{X}.
#' @usage d.spls.simulate(n=200,p=100,nondes=50,sigmaondes=0.05,sigmay=0.5,int.coef=1:5)
#' @param n a positive integer. \code{n} is the number of observations. Default value is \code{200}.
#' @param p a numeric vector of length \code{G} representing the number of variables. Default value is \code{100}.
#' @param nondes a numeric vector of length \code{G}. \code{nondes} is the number of Guassians in each mixture. Default value is \code{50}.
#' @param sigmaondes a numeric vector of length \code{G}. \code{sigmaondes} is the standard deviation of the
#' Gaussians for each group \eqn{g \in [1,G]}. Default value is \code{0.05}.
#' @param sigmay a real value. \code{sigmay} is the uncertainty on \code{y}. Default value is \code{0.5}.
#' @param int.coef a numeric vector of the coefficients of the linear combination in the construction of the response
#' vector \code{y}.
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
#' \deqn{y_i= \sigma_y \times V_i +\sum\limits_{g=1}^G \sum\limits_{k=1}^K \textrm{int.coeff}_k \times \textrm{sum}X^{g}_{ik}}
#' Where
#' \itemize{
#' \item \eqn{G} is the number of predictor sub matrices,
#' \item \eqn{i} is the index of the observation,
#' \item \eqn{V} is a normally distributed vector of 0 mean and unitary standard deviation,
#' \item \eqn{K} is the length of the vector \code{int.coeff},
#' \item \eqn{\textrm{sum}X^{g}} is a matrix of \eqn{n} rows and \eqn{K} columns.
#' The values of the column \eqn{k} are the sum of selected parts of each row of the sub matrix \eqn{X^g}. The columns of \eqn{X^g} are
#' seperated equally and each part is used for the \eqn{K} columns of \eqn{\textrm{sum}X^{g}}.
#' }
#'
#' @return A \code{list} of the following attributes
#' \item{X}{the concatunated predictors matrix.}
#' \item{y}{the response vector.}
#' \item{y0}{the reponse vector without noise \code{sigmay}.}
#' \item{sigmay}{the uncertainty on \code{y}.}
#' \item{sigmaondes}{the standard deviation of the Gaussians.}
#' \item{G}{the number of groups.}
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


d.spls.simulate<- function(n=200,p=100,nondes=50,sigmaondes=0.05,sigmay=0.5,int.coef=1:5)
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
  pif=round(seq(10,100,length.out = 2*length(int.coef))*p[1]/100)

  #Computing y0 as a sum of intervals of X
  sumX=matrix(0,n,length(int.coef))
  for (i in 1:length(int.coef))
  {
    sumX[,i]=apply(X[,pif[i]:pif[i+1]],1,function(u) sum(u))
  }
  y0=sumX%*%int.coef


  if (length(p)>1)
  {
    for (k in 2:length(p))
    {
      pif=round(seq(10,100,length.out = 2*length(int.coef))*p[k]/100)

      sumX=matrix(0,n,length(int.coef))
      for (i in 1:length(int.coef))
      {
        sumX[,i]=apply(X[,(sum(p[1:(k-1)])+pif[i]):(sum(p[1:(k-1)])+pif[i+1])],1,function(u) sum(u))
      }
      y0=y0+sumX%*%int.coef

    }
  }
  #adding noise to y0
  y=y0+sigmay*rnorm(n)
  y=as.vector(y)
  G=length(p)

  return(list(X=X,y=y,y0=y0,sigmay=sigmay,sigmaondes=sigmaondes,G=G))
}
