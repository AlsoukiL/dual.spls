#' Dual Sparse Partial Least Squares (Dual-SPLS) regression for the group lasso norm C
#' @keywords internal
#' @description
#' The function \code{d.spls.GLC} performs dimentional reduction as in PLS methodology combined to variable selection using the
#' Dual-SPLS algorithm with the norm \eqn{\Omega(w)=\|w_g\|_2+ \lambda_g \|w_g\|_1} for combined data where
#' \eqn{\Omega(w)=\sum\limits_{g=1}^G \alpha_g \Omega_g(w)=1; \sum\limits_{g=1}^G \alpha_g=1 \textrm{ and  } \forall g_1 \textrm{,  } g_2 \in \{1,\dots,G\}, \|w_{g_1} \|_2=\|w_{g_2} \|_2.}. Where \code{G} is the number of groups.
#' Dual-SPLS for the group lasso norms has been designed to confront the situations where the predictors
#' variables can be divided in distinct meaningful groups. Each group is constrained by an independent
#' threshold as in the dual sparse lasso methodology,
#' that is each \eqn{w_g} will be colinear to a vector \eqn{z_{\nu_g}} built from the coordinate of \eqn{z}
#' and constrained by the threshold \eqn{\nu_g}. Norm C applies the lasso norm for each group individually while constraining the overall norm.
#' @param X a numeric matrix of predictors values of dimension \code{(n,p)}. Each row represents one observation and each column one predictor variable.
#' @param y a numeric vector or a one column matrix of responses. It represents the response variable for each converstation.
#' @param ncp a positive integer. \code{ncp} is the number of Dual-SPLS components.
#' @param ppnu a positive real value or a vector of length the number of groups, in \eqn{[0,1]}.
#' \code{ppnu} is the desired proportion of variables to shrink to zero for each component and for each group.
#' @param indG a numeric vector of group index for each observation.
#' @param verbose a boolean value indicating whether or not to diplay the iterations steps.
#' @return A \code{list} of the following attributes
#' \item{Xmean}{the mean vector of the predictors matrix \code{X}.}
#' \item{scores}{the matrix of dimension \code{(n,ncp)} where \code{n} is the number of observations.The \code{scores} represents
#' the observations in the new component basis computed by the compression step
#' of the Dual-SPLS.}
#' \item{loadings}{the matrix of dimension \code{(p,ncp)} that represents the Dual-SPLS components.}
#' \item{Bhat}{the matrix of dimension \code{(p,ncp)} that regroups the regression coefficients for each component.}
#' \item{intercept}{the vector of length \code{ncp} representing the intercept values for each component.}
#' \item{fitted.values}{the matrix of dimension \code{(n,ncp)} that represents the predicted values of \code{y}}
#' \item{residuals}{the matrix of dimension \code{(n,ncp)} that represents the residuals corresponding
#'  to the difference between the responses and the fitted values.}
#' \item{lambda}{the matrix of dimension \code{(G,ncp)} collecting the parameters of sparsity \eqn{\lambda_g} used to fit the model at each iteration and for each group.}
#' \item{alpha}{the matrix of dimension \code{(G,ncp)} collecting the constraint parameters \eqn{\alpha_g}  used to fit the model at each iteration and for each group.}
#' \item{zerovar}{the matrix of dimension \code{(G,ncp)} representing the number of variables shrinked to zero per component and per group.}
#' @author Louna Alsouki François Wahl
#' @seealso [dual.spls::d.spls.GLA()],[dual.spls::d.spls.GLB()],[dual.spls::d.spls.GL()],`browseVignettes("dual.spls")`
#'
d.spls.GLC<- function(X,y,ncp,ppnu,indG,verbose=FALSE)
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
  listealpha=matrix(0,GG,ncp)


  nu=array(0,GG) #Initialising nu for each group
 # lambda=array(0,GG) #Initialising lambda for each group
#  alpha=array(0,GG) #Initialising alpha for each group
  Znu=array(0,p) #Initialising Znu for each group
  w=array(0,p) #Initialising w for each group
  norm2Znu=array(0,GG) #Initialising norm2 of Znu for each group
  norm1Znu=array(0,GG) #Initialising norm1 of Znu for each group

  ###################################
  # Dua-SPLS
  ###################################
  sam=10
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
      iz=which.min(abs(Zsp-ppnu[ig]))
      ###########
      nu[ig]=Zs[iz] #
      ###########

      # finding lambda, mu, given nu
      Znu[ind]=sapply(Z[ind],function(u) sign(u)*max(abs(u)-nu[ig],0))
      #Znu2=d.spls.norm2(Znu)
      #Znu1=d.spls.norm1(Znu)

      ##########Norm 1 of Znu(g)#############
      norm1Znu[ig]=d.spls.norm1(Znu[ind])
      ##########Norm 2 of Znu(g)#############
      norm2Znu[ig]=d.spls.norm2(Znu[ind])

    }
    #######################
    mu=sum(norm2Znu)

    #######################

    alpha=norm2Znu/mu
    lambda=nu/(mu*alpha)

    #Computing max of each norm2 of wg
    max_norm2w=array(0,GG) #Initialising maximum value of norm2 of w for each group

    for ( igg in 1:GG)
    {
      ######################
      max_norm2w[igg]=1/(alpha[igg]*(1+(nu[igg]*norm1Znu[igg]/(mu*alpha[igg])^2)))
      ######################
    }

    #Sampling the possible values of wg
    sample_wg=matrix(0,sam,(GG-1))
    for ( igg in 1:(GG-1))
    {
      ######################
      sample_wg[,igg]=seq(0,max_norm2w[igg],length.out = sam)
      ######################
    }

    #possible combination
    temp=data.frame(sample_wg[,1:(GG-1)])
    comb=expand.grid(temp)
    comb=unique(comb)
    comb=as.matrix(comb)
    combnum=dim(comb)[1]
    comb=cbind(comb,rep(0,combnum))
    denom=alpha[(GG)]*(1+( nu[(GG)] * norm1Znu[(GG)]  /  (mu*alpha[(GG)])^2   ))

    for (j in 1:combnum)
    {
      numg=rep(0,(GG-1))
      for ( igg in 1:(GG-1))
      {
        ######################
        numg[igg]=alpha[igg]*comb[j,igg]*(1+( nu[igg] * norm1Znu[igg]  /  (mu*alpha[igg])^2   ))
      }

      num=1-sum(numg)
      comb[j,GG]=num/denom
    }
    #verification
    # a=0
    # norm2wg=comb[14,]
    # for ( igg in 1:(GG))
    # {
    #   a=a+alpha[igg]*norm2wg[igg]*(1+( nu[igg] * norm1Znu[igg]  /  (mu*alpha[igg])^2   ))
    # }
    comb=comb[comb[,GG]>=0,]
    combnum=dim(comb)[1]
    RMSE=rep(0,combnum)
    tempw=matrix(0,p,combnum)
    tempt=matrix(0,n,combnum)
    tempBhat=matrix(0,p,combnum)
    tempintercept=rep(0,combnum)
    for (j in 1:combnum)
    {

      for ( igg in 1:(GG-1))
      {
        #Index of the group
        ind=which(indG==igg)
        # calculating w,t at the optimum
        w[ind]=(comb[j,igg]/(mu*alpha[igg]))*Znu[ind]
      }
      ind=which(indG==GG)
      w[ind]=(comb[j,GG]/(mu*alpha[GG]))*Znu[ind]

      #Finding T
      t=Xdef%*%w
      t=t/d.spls.norm2(t)

      WW[,ic]=w
      TT[,ic]=t

      #Coefficient vectors
      R=t(TT[,1:ic,drop=FALSE])%*%Xc%*%WW[,1:ic,drop=FALSE]
      R[row(R)>col(R)]<-0 # inserted for numerical stability

      L=backsolve(R,diag(ic))
      Bhat[,ic]=WW[,1:ic,drop=FALSE]%*%(L%*%(t(TT[,1:ic,drop=FALSE])%*%yc))

      tempintercept[ic] = ym - Xm %*% Bhat[,ic]

      #Predictions
      YY[,ic]=X %*% Bhat[,ic] + tempintercept[ic]
      RES[,ic]=y-YY[,ic]

      tempw[,j]=w
      tempt[,j]=t
      tempBhat[,j]=Bhat[,ic]
      RMSE[j]=sum(RES[,ic]^2)/n

    }
    indwmax=which.min(RMSE)
    w=tempw[,indwmax]
    t=tempt[,indwmax]

    WW[,ic]=w
    TT[,ic]=t

    Bhat[,ic]=tempBhat[,indwmax]

    listelambda[,ic]=lambda
    listealpha[,ic]=alpha
    intercept[ic] = ym - Xm %*% Bhat[,ic]

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
    #Deflation
    Xdef=Xdef-t%*%t(t)%*%Xdef


    # results iteration
    if (verbose){
      cat('Dual PLS ic=',ic,'lambda=',lambda,
          'mu=',mu,'nu=',nu,
          'nbzeros=',zerovar[,ic], '\n')
    }
  }

  return(list(Xmean=Xm,scores=TT,loadings=WW,Bhat=Bhat,intercept=intercept,
              fitted.values=YY,residuals=RES,
              lambda=listelambda,alpha=listealpha,zerovar=zerovar))
}
