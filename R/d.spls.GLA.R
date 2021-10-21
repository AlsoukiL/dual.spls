#' Dual Sparse Partial Least Squares (Dual-SPLS) regression for the norm \eqn{\Omega(w)=\|w\|_2+\sum_{g} \lambda_g\|w_g\|_1}
#' @keywords internal
#' @description
#' The function \code{d.spls.GLA} performs dimentional reduction combined to variable selection using the
#' Dual-SPLS algorithm with the norm \eqn{\Omega(w)=\|w\|_2+\sum_{G} \lambda_g\|w_g\|_1} for combined data.
#' @param X a numeric matrix of predictors values. Each row represents an observation and each column a predictor variable.
#' @param y a numeric vector or a one column matrix of responses. It represents the response variable for each converstation.
#' @param ncp a positive integer. \code{ncp} is the number of Dual-SPLS components.
#' @param pctnu a positive real value or a vector of length the number of groups, in \code{[0,1]}.
#' \code{pctnu} is the desired proportion of variables to shrink to zero for each component and for each group.
#' @param indG a numeric vector of group index for each observation.
#' @param verbose a boolean value indicating whether or not to diplay the iterations steps.
#' @return A \code{list} of the following attributes
#' \itemize{
#' \item Xmean the mean vector of the predictors matrix \code{X}.
#' \item scores a matrix of \code{ncp} columns representing the Dual-SPLS components and the same number of rows
#' as \code{X} representing the observations in the new component basis computed by the compression step
#' of the Dual-SPLS.
#' \item loadings the loadings matrix.
#' \item Bhat the regression coefficients matrix for each component.
#' \item intercept the intercept value for each component.
#' \item fitted.values the matrix of predicted values of \code{y}
#' \item residuals the matrix of residuals corresponding to the difference between the response and the fitted values.
#' \item lambda the matrix of the sparse hyper-parameter used to fit the model at each iteration and for each group.
#' \item zerovar the matrix of the number of variables shrinked to zero per component and per group.
#' }
#' @seealso [dual.spls::d.spls.GLB()],[dual.spls::d.spls.GLC()],[dual.spls::d.spls.GL()],`browseVignettes("dual.spls")`
#'
#'
#'
d.spls.GLA<- function(X,y,ncp,pctnu,indG,verbose=FALSE)
{

  ###################################
  # Dimensions
  ###################################
  n=length(y) #Number of observations
  p=dim(X)[2] #Number of variables

  ###################################
  # Centering Data
  ###################################
  Xm = apply(X, 2, mean) #Mean of X
  Xc=X - rep(1,n) %*% t(Xm) #Centering predictor matrix

  ym=mean(y) #Mean of y
  yc=y-ym #Centering response vector

  ###################################
  # Initialisation
  ###################################
  GG=max(indG) #Number of groups
  PP=array(0,GG) #Initialising the vector of dimension of each group


  for (ig in 1:GG)
  {
    PP[ig]=sum(indG==ig)
  }


  WW=matrix(0,p,ncp) #Initialising W, the matrix of loadings
  TT=matrix(0,n,ncp) #Initialising T, the matrix of scores
  Bhat=matrix(0,p,ncp) #Initialising the matrix of coefficients
  YY=matrix(0,n,ncp) #Initialising the matrix of coefficients
  RES=matrix(0,n,ncp) #Initialising the matrix of coefficients
  intercept=rep(0,ncp)
  zerovar=matrix(0,GG,ncp)
  listelambda=matrix(0,GG,ncp)


  nu=array(0,GG) #Initialising nu for each group
  lambda=array(0,GG) #Initialising lambda for each group
  Znu=array(0,p) #Initialising Znu for each group
  w=array(0,p) #Initialising w for each group
  norm2Znu=array(0,GG) #Initialising norm2 of Znu for each group
  norm1Znu=array(0,GG) #Initialising norm1 of Znu for each group

  ###################################
  # Dual-SPLS
  ###################################

  #Each step ic in -for loop- determine the icth column of each W, T and Bhat
  Xdef=Xc #Initialising X for Deflation Step
  for (ic in 1:ncp)
  {

    Z=t(Xdef)%*%yc #For cov(t(X)y,w)=0, w must be colinear to t(X)y ==> Z=t(X)y
    Z=as.vector(Z)

    for( ig in 1:GG)
    {
      #Index of the group
      ind=which(indG==ig)

      #Optimizing nu(g)
      Zs=sort(abs(Z[ind]))
      d=length(Zs)
      Zsp=(1:d)/d
      iz=which.min(abs(Zsp-pctnu[ig]))
      ###########
      nu[ig]=Zs[iz] #
      ###########

      # finding lambda, mu, given nu
      Znu[ind]=sapply(Z[ind],function(u) sign(u)*max(abs(u)-nu[ig],0))
      #Znu2=norm2(Znu)
      #Znu1=norm1(Znu)

      ##########Norme 1 de Znu(g)#############
      norm1Znu[ig]=norm1(Znu[ind])
      ##########Norme 2 de Znu(g)#############
      norm2Znu[ig]=norm2(Znu[ind])

    }
    #######################
    mu=norm2(Znu)
    mu2=mu^2

    #######################

    for ( igg in 1:GG)
    {
      ######################
      lambda[igg]=nu[igg]/mu #
      ######################
      #Index of the group
      ind=which(indG==igg)

      # calculating w,t at the optimum
      w[ind]=(mu/(mu2+ t(nu)%*%norm1Znu))%*%Znu[ind]
    }

    #Finding T
    t=Xdef%*%w
    t=t/norm2(t)

    WW[,ic]=w
    TT[,ic]=t

    #Deflation
    Xdef=Xdef-t%*%t(t)%*%Xdef

    #Coefficient vectors
    R=t(TT[,1:ic,drop=FALSE])%*%Xc%*%WW[,1:ic,drop=FALSE]
    R[row(R)>col(R)]<-0 # inserted for numerical stability

    L=backsolve(R,diag(ic))
    Bhat[,ic]=WW[,1:ic,drop=FALSE]%*%(L%*%(t(TT[,1:ic,drop=FALSE])%*%yc))

    listelambda[,ic]=lambda
    intercept[ic] = ym - Xm %*% Bhat[,ic]

    #Number of zero variables in each group
    zerovar[1,ic]=sum(Bhat[(1:PP[1]),ic]==0)
    if (GG>1){
      for (iggg in 1:(GG-1))
      {
        zerovar[iggg+1,ic]=sum(Bhat[((PP[iggg]+1):(PP[iggg]+PP[iggg+1])),ic]==0)
      }
    }
    #Predictions
    YY[,ic]=X %*% Bhat[,ic] + intercept[ic]
    RES[,ic]=y-YY[,ic]

    # results iteration
    if (verbose){
      cat('Dual PLS ic=',ic,'lambda=',lambda,
          'mu=',mu,'nu=',nu,
          'nbzeros=',zerovar[,ic], '\n')
    }
  }

  return(list(Xmean=Xm,scores=TT,loadings=WW,Bhat=Bhat,intercept=intercept,
              fitted.values=YY,residuals=RES,
              lambda=listelambda,zerovar=zerovar))
}
