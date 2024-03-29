#' Dual Sparse Partial Least Squares (Dual-SPLS) regression for the group lasso norm A
#' @keywords internal
#' @description
#' The function \code{d.spls.GLA} performs dimensional reduction as in PLS methodology combined to variable selection using the
#' Dual-SPLS algorithm with the norm \deqn{\Omega_g(w)=\|w_g\|_2+ \lambda_g \|w_g\|_1} for combined data where
#' \eqn{\Omega(w)=\sum_{g=1}{^G} \alpha_g \Omega_g(w)=1};
#' \eqn{\sum_{g=1}^G \alpha_g=1} and \code{G} is the number of groups.
#' Dual-SPLS for the group lasso norms has been designed to confront the situations where the predictors
#' variables can be divided in distinct meaningful groups. Each group is constrained by an independent
#' threshold as in the dual sparse lasso methodology,
#' that is each \eqn{w_g} will be collinear to a vector \eqn{z.\nu_g} built from the coordinate of \eqn{z}
#' and constrained by the threshold \eqn{\nu_g}. Norm A i a generalized group lasso-like norm that applies the lasso norm for each group individually while constraining the overall norm. Moreover,
#' the Euclidian norm of each \eqn{w_g} is computed while minimizing the root mean squares error of prediction.
#' @param X a numeric matrix of predictors values of dimension \code{(n,p)}. Each row represents one observation and each column one predictor variable.
#' @param y a numeric vector or a one column matrix of responses. It represents the response variable for each observation.
#' @param ncp a positive integer. \code{ncp} is the number of Dual-SPLS components.
#' @param ppnu a positive real value or a vector of length the number of groups, in \eqn{[0,1]}.
#' \code{ppnu} is the desired proportion of variables to shrink to zero for each component and for each group.
#' @param indG a numeric vector of group index for each observation.
#' @param verbose a Boolean value indicating whether or not to display the iterations steps.
#' @return A \code{list} of the following attributes
#' \item{Xmean}{the mean vector of the predictors matrix \code{X}.}
#' \item{scores}{the matrix of dimension \code{(n,ncp)} where \code{n} is the number of observations. The \code{scores} represents
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
#' \item{zerovar}{the matrix of dimension \code{(G,ncp)} representing the number of variables shrank to zero per component and per group.}
#' \item{PP}{the vector of length \code{G} specifying the number of variables in each group.}
#' \item{ind_diff0}{the list of \code{ncp} elements representing the index of the none null regression coefficients elements.}
#' \item{type}{a character specifying the Dual-SPLS norm used. In this case it is \code{GLA}. }
#' @author Louna Alsouki François Wahl
#' @seealso [dual.spls::d.spls.GLA],[dual.spls::d.spls.GLB],[dual.spls::d.spls.GL]
#'
d.spls.GLA<- function(X,y,ncp,ppnu,indG,verbose=FALSE)
{

  ###################################
  # Dimensions
  ###################################
  n=length(y) # number of observations
  p=dim(X)[2] # number of variables

  ###################################
  # Centering Data
  ###################################
  Xm = apply(X, 2, mean) # mean of X
  Xc=X - rep(1,n) %*% t(Xm) # centering predictor matrix

  ym=mean(y) # mean of y
  yc=y-ym # centering response vector

  ###################################
  # Initialisation
  ###################################
  nG=max(indG) # number of groups
  PP=sapply(1:nG, function(u) sum(indG==u) )

  WW=matrix(0,p,ncp) # initializing WW, the matrix of loadings
  TT=matrix(0,n,ncp) # initializing TT, the matrix of scores
  Bhat=matrix(0,p,ncp) # initializing Bhat, the matrix of coefficients
  YY=matrix(0,n,ncp) # initializing YY, the matrix of coefficients
  RES=matrix(0,n,ncp) # initializing RES, the matrix of coefficients
  intercept=rep(0,ncp) # initializing intercept, the vector of intercepts
  zerovar=matrix(0,nG,ncp) # initializing zerovar, the matrix of final number of zeros coefficients for each component and for each group
  listelambda=matrix(0,nG,ncp) # initializing listelambda, the matrix of values of lambda for each group
  listealpha=matrix(0,nG,ncp) # initializing listealpha, the matrix of values of alpha for each group
  ind.diff0=vector(mode = "list", length = ncp) # initializing ind0, the list of the index of the none zero coefficients
  names(ind.diff0)=paste0("in.diff0_", 1:ncp)

  nu=array(0,nG) # initializing nu for each group
  Znu=array(0,p) # initializing Znu for each group
  w=array(0,p) # initializing w for each group
  norm2Znu=array(0,nG) # initializing norm2 of Znu for each group
  norm1Znu=array(0,nG) # initializing norm1 of Znu for each group

  ###################################
  # Dua-SPLS
  ###################################
  # each step ic in -for loop- determine the icth column or element of each element initialized
  Xdef=Xc # initialising X for Deflation Step
  for (ic in 1:ncp)
  {
    Z=t(Xdef)%*%yc
    Z=as.vector(Z)

    for( ig in 1:nG)
    {
      # index of the group
      ind=which(indG==ig)

      # optimizing nu(g)
      Zs=sort(abs(Z[ind]))
      d=length(Zs)
      Zsp=(1:d)/d
      iz=which.min(abs(Zsp-ppnu[ig]))
      ###########
      nu[ig]=Zs[iz] #
      ###########

      # finding mu, given nu
      Znu[ind]=sapply(Z[ind],function(u) sign(u)*max(abs(u)-nu[ig],0))

      ##########Norm 1 of Znu(g)#############
      norm1Znu[ig]=d.spls.norm1(Znu[ind])
      ##########Norm 2 of Znu(g)#############
      norm2Znu[ig]=d.spls.norm2(Znu[ind])
    }
    #######################
    mu=sum(norm2Znu)
    #######################

    # finding alpha and lambda, given nu
    alpha=norm2Znu/mu
    lambda=nu/(mu*alpha)

    # computing max of each norm2 of wg
    max_norm2w=sapply(1:nG,function(u) 1/alpha[u]/(1+(nu[u]*norm1Znu[u]/(mu*alpha[u])^2)))

    # sampling the possible values of wg
    sample_wg=matrix(0,10,(nG-1))
    for ( ig in 1:(nG-1))
    {
      ######################
      sample_wg[,ig]=seq(0,max_norm2w[ig],length.out = 10)
      ######################
    }
    # all the possible combinations
    # due to a linear constraint, the last column is deducted from the others
    temp=data.frame(sample_wg[,1:(nG-1)])
    comb=expand.grid(temp)
    comb=unique(comb)
    comb=as.matrix(comb)
    ncomb=dim(comb)[1]
    comb=cbind(comb,rep(0,ncomb))
    denom=alpha[(nG)]*(1+( nu[(nG)] * norm1Znu[(nG)]  /  (mu*alpha[(nG)])^2   ))
    for (icomb in 1:ncomb)
    {
      numg=sapply(1:(nG-1), function(u){
        alpha[u]*comb[icomb,u]*(1+( nu[u]*norm1Znu[u]/(mu*alpha[u])^2)) } )
      num=1-sum(numg)
      comb[icomb,nG]=num/denom
    }
    # suppress inadequate rows (for which comb[,nG] < 0)
    comb=comb[comb[,nG]>=0,]
    ncomb=dim(comb)[1]

    # initializing RMSE
    RMSE=rep(0,ncomb)

    # initializing RMSE
    tempw=matrix(0,p,ncomb)

    # computing w for each combination
    for (icomb in 1:ncomb)
    {

      for ( ig in 1:nG )
      {
        #Index of the group
        ind=which(indG==ig)
        # calculating w,t at the optimum
        w[ind]=(comb[icomb,ig]/(mu*alpha[ig]))*Znu[ind]
      }

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

      intercept[ic] = ym - Xm %*% Bhat[,ic]

      #Predictions
      YY[,ic]=X %*% Bhat[,ic] + intercept[ic]
      RES[,ic]=y-YY[,ic]

      tempw[,icomb]=w
      RMSE[icomb]=sum(RES[,ic]^2)/n

    }

    # choosing the optimal w
    indwmax=which.min(RMSE)
    w=tempw[,indwmax]
    WW[,ic]=w

    # finding T
    t=Xdef%*%w
    t=t/d.spls.norm2(t)
    TT[,ic]=t

    #deflation
    Xdef=Xdef-t%*%t(t)%*%Xdef

    # coefficient vectors
    R=t(TT[,1:ic,drop=FALSE])%*%Xc%*%WW[,1:ic,drop=FALSE]
    R[row(R)>col(R)]<-0 # inserted for numerical stability

    L=backsolve(R,diag(ic))
    Bhat[,ic]=WW[,1:ic,drop=FALSE]%*%(L%*%(t(TT[,1:ic,drop=FALSE])%*%yc))

    # lambda
    listelambda[,ic]=lambda

    # alpha
    listealpha[,ic]=alpha

    # intercept
    intercept[ic] = ym - Xm %*% Bhat[,ic]

    # zerovar
    zerovar[,ic]=sapply(1:nG, function(u) {
      indu=which(indG==u)
      sum(Bhat[indu,ic]==0)})

    # ind.diff0
    name=paste("in.diff0_",ic,sep="")
    ind.diff0[name]=list(which(Bhat[,ic]!=0))

    # predictions
    YY[,ic]=X %*% Bhat[,ic] + intercept[ic]

    # residuals
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
              lambda=listelambda,alpha=listealpha,zerovar=zerovar,PP=PP,ind.diff0=ind.diff0,type="GLA"))
}
