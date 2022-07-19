#' Print function for Dual-SPLS models
#' @description
#' The function \code{d.spls.print} pisplays the values of a Dual-SPLS regression parameters of sparsity and the number of variables
#' shrinked to zero for a specified number of components.
#' @usage d.spls.print(mod.dspls,ncomp)
#' @param mod.dspls is a fitted Dual-SPLS object.
#' @param ncomp a positive integer of the number of Dual-SPLS components.
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
#' ncomplasso <- d.spls.cv(X=X,Y=y,ncomp=10,dspls="lasso",ppnu=0.9,nrepcv=20,pctcv=75)
#' mod.dspls.lasso <- d.spls.lasso(X=X,y=y,ncp=ncomplasso,ppnu=0.9,verbose=TRUE)
#'
#' d.spls.print(mod.dspls.lasso,ncomplasso)
#' @export

d.spls.print<-function (mod.dspls,ncomp)
{
  Xmean=mod.dspls$Xmean
  lambda=mod.dspls$lambda[ncomp]
  zerovar=mod.dspls$zerovar[ncomp]
  type=mod.dspls$type
  p=length(Xmean)
  cat( "\nDual Sparse Partial Least Squares Regression for the ", type, " norm","\n" )
  cat( "--------------------------------------------\n\n")
  if (type=="GLB" || type=="GLC") {
    lambda=mod.dspls$lambda[,ncomp]
    zerovar=mod.dspls$zerovar[,ncomp]
    alpha=mod.dspls$alpha[,ncomp]
    PP=mod.dspls$PP
    for(i in 1:length(PP)){
      cat( paste("For group ",i," : \n"))
      cat( paste("Parameters: lambda =",lambda[i],", alpha = ",alpha[i],", ncomp = ",ncomp,"\n") )
      cat( paste("Dual-SPLS selected: ",zerovar[i]," variables among ",PP[i]," variables from group ",i," \n\n") )
    }
  }
  else{
    if (type=="GLA") {
      lambda=mod.dspls$lambda[,ncomp]
      zerovar=mod.dspls$zerovar[,ncomp]
      alpha=mod.dspls$alpha[,ncomp]
      PP=mod.dspls$PP
      for(i in 1:length(PP)){
        cat( paste("For group ",i," : \n"))
        cat( paste("Parameters: lambda =",lambda[i],", ncomp = ",ncomp,"\n") )
        cat( paste("Dual-SPLS selected: ",zerovar[i]," variables among ",PP[i]," variables from group ",i," \n\n") )
      }
    }
    else
    {
      cat( paste("Parameters: lambda = ",lambda,", ncomp = ",ncomp,"\n") )
      cat( paste("Dual-SPLS selected: ",zerovar," variables among ",p," variables"," \n") )

    }
  }
}
