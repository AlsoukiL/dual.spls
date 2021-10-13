#' Dual Sparse Partial Least Squares (Dual-SPLS) regression for the group lasso norms
#' @description
#' The function \code{d.spls.GL} performs dimentional reduction combined to variable selection using the
#' Dual-SPLS algorithm with the group lasso norms
#' \itemize{
#' \item Norm A: \eqn{\Omega(w)=\|w\|_2+\sum_{g} \lambda_g\|w_g\|_1}
#' \item Norm B: \eqn{\Omega(w)=\sum_{g} \alpha_g \|w \|_2+\sum_{g} \lambda_g \|w_g \|_1} for
#' \eqn{\sum_{g} \alpha_g=1; \Omega(w_g)=\gamma_g ;\sum_{g} \gamma_g=1}
#' \item Norm C: \eqn{\Omega(w)=\|w_g\|_2+ \lambda_g \|w_g\|_1} for
#' \eqn{\Omega(w)=\sum_{g} \alpha_g \Omega_g(w)=1; \sum_{g} \alpha_g=1}
#' }
#' @usage
#' \code{d.spls.GL(X,y,ncp,pctnu,G,gamma=NULL,norm="A",verbose=FALSE) }
#' @param X a numeric matrix of predictors values. Each row represents an observation and each column a predictor variable.
#' @param y a numeric vector or a one column matrix of responses. It represents the response variable for each converstation.
#' @param ncp a positive integer. \code{ncp} is the number of Dual-SPLS components.
#' @param pctnu a positive real value or a vector of length the number of groups, in \code{[0,1]}.
#' \code{pctnu} is the desired proportion of variables to shrink to zero for each component and for each group.
#' @param G a numeric vector of group index for each variable.
#' @param gamma a numeric vector of the norm \eqn{\Omega} of each \eqn{w_g} in case \code{norm="A"}. Default value is NULL.
#' @param norm a character specifying the norm chosen between A, B and C. Default value is A.
#' @param verbose a boolean value indicating whether or not to diplay the iterations steps.
#' @details
#' This procdure computes latent sparse components that are used in a regression model.
#' The optimization problem of the Dual-SPLS regression for the group lasso norms is
#' similar to the Dual norm defintion of the norm \eqn{\Omega(w)}. Indeed, we are searching for \code{w} that goes
#' with
#' \deqn{\Omega^*(z)=max_w(z^Tw) \text{ s.t. } \Omega(w)=1}
#' The solution of this problem for the norm A is
#' \deqn{\dfrac{w_g}{\|w\|_2}=\dfrac{1}{\mu} \delta_g (|z_g|-\nu_g)_+ \text{ for each group} g}
#'
#' The solution of this problem for the norm B is
#' \deqn{\dfrac{w_g}{\|w\|_2}=\dfrac{1}{\mu \alpha_g} \delta_g (|z_g|-\nu_g)_+ \text{ for each group} g}
#'
#' The solution of this problem for the norm C is
#' \deqn{\dfrac{w_g}{\|w_g\|_2}=\dfrac{1}{\mu \alpha_g} \delta_g (|z_g|-\nu_g)_+ \text{ for each group} g}
#'
#'
#' Where
#' \itemize{
#' \item \eqn{\delta_g} is the vector of signs of \code{z_g} and \code{w_g}
#' \item \eqn{\alpha_g} is the hyper-parameter of the procedure.
#' \item \eqn{\mu} is a parameter that guarentees the constraint of \eqn{\Omega(w)=1}
#' \item \eqn{\nu_g} is the shrinkage parameter for each group \code{g}.
#' }
#'
#' The parameters \eqn{\nu_g} are computed adaptively according to the proportion of zero variables desired \code{pctnu} for each group.
#' @return A \code{list} of the following attributes
#' @param Xmean the mean vector of the predictors matrix \code{X}.
#' @param scores a matrix of \code{ncp} columns representing the Dual-SPLS components and the same number of rows
#' as \code{X} representing the observations in the new component basis computed by the compression step
#' of the Dual-SPLS.
#' @param loadings the loadings matrix.
#' @param Bhat the regression coefficients matrix for each component.
#' @param intercept the intercept value for each component.
#' @param fitted.values the matrix of predicted values of \code{y}
#' @param residuals the matrix of residuals corresponding to the difference between the response and the fitted values.
#' @param lambda the matrix of the first sparse hyper-parameter used to fit the model at each iteration and for each group.
#' @param alpha the matrix of the second sparse hyper-parameter used to fit the model at each iteration and for each group
#' when the norm chosen is B or C.
#' @param zerovar the matrix of the number of variables shrinked to zero per component and per group.
#' @seealso [dual.spls::d.spls.GLA()],[dual.spls::d.spls.GLB()],[dual.spls::d.spls.GLC()],`browseVignettes("dual.spls")`
#'
#'
#' @examples
#' ### load dual.spls library
#' library(dual.spls)
#'
#' @export


d.spls.GL<- function(X,y,ncp,pctnu,G,gamma=NULL,norm="A",verbose=FALSE)
{
  if (norm=="A")
  {
    mod.dspls=DUALSPLSgrouplassoA(X=X,y=y,ncp=ncp,pctnu=pctnu,G=G,verbose=verbose)
  }
  if (norm=="B")
  {
    mod.dspls=DUALSPLSgrouplassoB(X=X,y=y,ncp=ncp,pctnu=pctnu,G=G,gamma=gamma,verbose=verbose)
  }

  if (norm=="C")
  {
    mod.dspls=DUALSPLSgrouplassoC(X=X,y=y,ncp=ncp,pctnu=pctnu,G=G,verbose=verbose)
  }
  return(mod.dspls)
}
