#' Computes the error of a cross validation iteration
#' @keywords internal
#' @description
#' The function \code{d.spls.errorcv} computes the error of a validation set according to a calibration set \code{cvcal} used
#' to fit the Dual-SPLS regression. This function is an internal function used in the cross validation procedure in order to determine
#' the best number of latent variables of any of the Dual-SPLS versions.
#' @param cvcal a numeric vector representing the index of the calibration set to be used in the fitting.
#' @param X a numeric matrix.
#' @param Y a numeric vector representing the response values.
#' @param ncomp a numeric vector of the number of latent numbers to use while computing the errors.
#' @param dspls the norm type of the Dual-SPLS regression applied. Default value is \code{lasso}. Options are \code{pls}, \code{LS},
#' \code{ridge}, \code{GLA}, \code{GLB} and \code{GLC}.
#' @param ppnu a positive real value, in \eqn{[0,1]}. \code{ppnu} is the desired
#' proportion of variables to shrink to zero for each component (see Dual-SPLS methodology).
#' @param nu2 a positive real value. \code{nu2} is a constraint parameter used in the ridge norm.
#' @param indG a numeric vector of group index for each observation. It is used in the cases of the group lasso norms.
#' @param gamma a numeric vector of the norm \eqn{\Omega} of each \eqn{w_g} in the case of \code{GLB} norm.
#' @return a numeric vector representing the errors for each fitted model
#' @author Louna Alsouki François Wahl
#' @seealso [dual.spls::d.spls.cv()],[dual.spls::d.spls.lasso()],`browseVignettes("dual.spls")`
#'
#'
d.spls.errorcv<-function (cvcal, X, Y, ncomp,dspls="lasso",ppnu=0.9,nu2,indG,gamma)
{

  if(dspls=="pls") mod.dpls=d.spls.pls(X[cvcal, ], Y[cvcal ], ncp = max(ncomp))
  if(dspls=="lasso") mod.dpls=d.spls.lasso(X[cvcal, ], Y[cvcal ], ncp = max(ncomp),ppnu)
  if(dspls=="LS") mod.dpls=d.spls.LS(X[cvcal, ], Y[cvcal ], ncp = max(ncomp),ppnu)
  if(dspls=="ridge") mod.dpls=d.spls.ridge(X[cvcal, ], Y[cvcal ], ncp = max(ncomp),ppnu,nu2)
  if(dspls=="GLA") mod.dpls=d.spls.GLA(X[cvcal, ], Y[cvcal ],ncp=max(ncomp),ppnu,indG)
  if(dspls=="GLB") mod.dpls=d.spls.GLB(X[cvcal, ], Y[cvcal ],ncp=max(ncomp),ppnu,indG,gamma)
  if(dspls=="GLC") mod.dpls=d.spls.GLC(X[cvcal, ], Y[cvcal ],ncp=max(ncomp),ppnu,indG)
  yhatval=d.spls.predict(mod.dpls,X[-cvcal, ],liste_cp=ncomp)
  errorcv=apply(yhatval,2,function(u) sum((u - Y[-cvcal ])^2))

  return(errorcv)
}
