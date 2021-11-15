#' @name dual.spls-package
#' @docType package
#' @title dual.spls package
#'
#' @description
#' Provides a series of functions that compute latent sparse components that are used in a regression model.
#' The optimization problem of the Dual-SPLS regression depends on the norm \eqn{\Omega(w)} chosen.
#' This procedure computes latent sparse components that are used in a regression model.
#' Indeed, we are searching for \code{w} that goes with
#' \deqn{\Omega^*(z)=\max\limits_{w} (z^Tw)  \textrm{ s.t. } \Omega(w)=1}
#' It also suggest a calibration and validation method based on a modified version of the Kennard
#'  and Stone Algorithm and a function that simulates data composed of Guassians mixtures.
#' @author Louna Alsouki François Wahl
#' @seealso [dual.spls::d.spls.lasso()], [dual.spls::d.spls.LS()], [dual.spls::d.spls.ridge()],
#' [dual.spls::d.spls.GL()], `browseVignettes("dual.spls")`
#'
#' @importFrom pdist pdist
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @import mathjaxr
"_PACKAGE"
