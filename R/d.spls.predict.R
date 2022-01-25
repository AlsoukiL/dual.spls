#' Makes predictions from a fitted Dual-SPLS model
#' @description
#' The function \code{d.spls} makes predictions from a fitted Dual-SPLS model.
#' @usage d.spls.predict(mod.dspls,X,liste_cp)
#' @param mod.dspls a fitted Dual-SPLS object.
#' @param X a numeric matrix of predictors values. Each row represents an observation and each column a predictor variable.
#' @param liste_cp a numeric vector of the components for which prediction is required.
#' @details
#' The coefficients computed in the Dual-SPLS object are used to predict the fitted response values of new matrix \code{X}.
#' Users can choose how many Dual-SPLS components should be used.
#' @return Vector or matrix of estimated responses.
#' @author François Wahl Louna Alsouki
#'
#' @examples
#' ### load dual.spls library
#' library(dual.spls)
#' ### parameters
#' n <- 100
#' p <- 50
#' nondes <- 20
#' sigmaondes <- 0.5
#' data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
#'
#' X <- data$X
#' y <- data$y
#'
#' pcal <- 70
#' ncells <- 3
#'
#' split <- d.spls.calval(X=X,pcal=pcal,y=y,ncells=ncells)
#'
#' indcal= split$indcal
#' indval= split$indval
#'
#' Xcal=X[indcal,]
#' Xval=X[indval,]
#' ycal=y[indcal]
#' yval=y[indval]
#'
#' #fitting the model
#' ncp=10
#' mod.dspls <- d.spls.lasso(X=Xcal,y=ycal,ncp=ncp,ppnu=0.9,verbose=TRUE)
#'
#' ycalhat=mod.dspls$fitted.values
#' rescal=mod.dspls$residuals
#' # predictions on validation
#' yvalhat=d.spls.predict(mod.dspls,Xval, liste_cp=1:ncp)
#'
#' #Computing RMSE error
#' RMSEcal.dspls=apply(rescal,2,function(u) sqrt(mean(u^2)) )
#' RMSEval.dspls=apply(yvalhat,2,function(u) sqrt(mean((yval-u)^2)) )
#'
#'
#'
#' @export

d.spls.predict<- function(mod.dspls,X,liste_cp)
{
  # dimensions
  n=dim(X)[1]

  ###################################
  # Centering Data
  ###################################
  Xm = mod.dspls$Xmean #Mean of X
  Xc=X - rep(1,n) %*% t(Xm) #Centering predictor matrix

  ncp=length(liste_cp)

  # predictions on X
  M1=matrix(mod.dspls$intercept[liste_cp],ncp,n)
  M1=t(M1)
  yhat=X %*% mod.dspls$Bhat[,liste_cp] + M1

  return(yhat)

}
