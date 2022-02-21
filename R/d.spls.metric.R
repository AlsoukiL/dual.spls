#' Computes predictions criterias
#' @description
#' This function computes evaluation metrics commonly used in modelling. It provides the values of the root mean square
#' error (RMSE), the mean absolute error (MAE) and the Rsquared.
#' @usage d.spls.metric(hat.y,real.y)
#' @param hat.y a numeric vector. It represents the fitted response variable for each observation using a Dual-SPLS method.
#' @param real.y a numeric vector. It represents the response variable for each observation.
#' @details
#' The Root Mean Square Error (RMSE) is the standard deviation of the residuals. It is computed as follows:
#' \deqn{RMSE=\sqrt{\sum\limits_{i=1}^n \frac{(y_{fi}-y_{ri})^2}{n}}.}{RMSE=\sqrt(1/n\sum (y.f-y.r)^2).}
#'
#' The Mean Absolute Error measures the average magnitude of the errors in a set of predictions,
#' without considering their direction. It is computed as follows:
#' \deqn{MAE=\frac{1}{n}\sum\limits_{i=1}^n |y_{fi}-y_{ri}|.}{MAE=1/n \sum |y_f-y_r|.}
#'
#' The Rsquared represents the proportion of the variance for a dependent
#' variable that's explained by an independent variable in a regresison. It is computed as follows:
#'
#'\deqn{R2=\frac{\sum\limits_{i=1}^n (y_{fi}-\bar{y})^2}{\sum\limits_{i=1}^n (y_{ri}-\bar{y})^2}.}{R2=\sum_{i=1,n} (y_f-bar.y)^2/ \sum_{i=1,n} (y_r-bar.y)^2.}
#'Where \eqn{\bar{y}=\frac{1}{n}\sum\limits_{i=1}^n y_{fi}}{bar.y=1/n\sum y_f}. Note that \eqn{y_f}{y.f} are the fitted values and \eqn{y_r}{y.r} are the real ones.
#'
#' @return A \code{list} of the following attributes
#' \item{RMSE}{the vector of the root mean square error values for each component.}
#' \item{MAE}{the vector of the mean absolute error values for each component.}
#' \item{Rsquared}{the vector of the Rsquared values for each component.}
#'
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
#' ncomplasso <- d.spls.cv(X=X,Y=y,ncomp=10,dspls="lasso",ppnu=0.9,nrepcv=20,pctcv=75)
#' mod.dspls.lasso <- d.spls.lasso(X=X,y=y,ncp=ncomplasso,ppnu=0.9,verbose=TRUE)
#'
#' predmetric= d.spls.metric(mod.dspls.lasso$fitted.values,y)
#'
#' #Error plots
#' plot(1:ncomplasso,predmetric$RMSE,
#' main=("Root mean squares error values"),xlab='Number of components',ylab='Errors',col='blue',pch=19)
#' lines(1:ncomplasso,predmetric$RMSE,col='blue')
#' points(1:ncomplasso,predmetric$MAE,col='red',pch=19)
#' lines(1:ncomplasso,predmetric$MAE,col='red')
#' points(1:ncomplasso,predmetric$R2,col='green',pch=19)
#' lines(1:ncomplasso,predmetric$R2,col='green')
#' legend("topright", legend = c("RMSE", "MAE", "R2"), bty = "n",
#'  cex = 0.8, col = c("blue", "red","green"), lty = c(1,1,1))
#' @export

d.spls.metric<-function (hat.y,real.y)
{
  #residuals
  res=real.y-hat.y

  #fitted values
  fitted.y=hat.y

  #RMSE
  RMSE=apply(res,2,function(u) sqrt(mean(u^2)) )

  #MAE
  MAE=apply(res,2,function(u) mean(abs(u)) )

  #R2
  ybar=mean(real.y) # mean of y
  SCT=sum((real.y-ybar)^2) # total sum of squares
  p=dim(res)[2] # number of components
  SSres=rep(0,p) # regression sum of squares
  for (i in 1:p) {SSres[i]=sum((res[,i])^2)}
  R2=1-SSres/SCT

  return(list(RMSE=RMSE,MAE=MAE,Rsquared=R2))
}
