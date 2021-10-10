#' Splits data into calibration and validation sets using the modified Kennard and Stone method
#' @description
#' The function \code{modified.KS} divides the data \code{X} into a calibration and a validation set using
#' the modified Kennard and Stone startegy.
#' @usage
#' \code{calval(X,pcal,Datatype=NULL,y=NULL,ncells=10) }
#' @param X a numeric matrix of predictors values.
#' @param pcal a positive integer. \code{p} is the number of calibration samples to be selected.
#' @param Datatype A vecor of index specifying each observation belonging to wich group index.
#' Default value is NULL, meaning the function will use the function \code{type} to compute the vector for \code{ncells}.
#' If NULL, parameter \code{y} should be specified.
#' @param y a numeric vector of responses. Default value is NULL, meaning as long as \code{Datatype} is specified,
#' \code{y} is not necessary.
#' @param ncells a positive integer. \code{ncells} is the number of groups dividing the observations.
#' @details
#' The modified Kennard and Stone algorithm allows to select samples using the classical Kennard and stone on
#' each group of observations one by one. It starts by selecting the point that is the furthest away from the centroid.
#' This point is assigned as the calibration set and is removed from the list of candidates. Then, it adentify to which
#' group belongs this first observation and consider the group that comes after.
#' #' @return A \code{list} of the following attributes
#' @param X the preditors matrix. It represent the data to split into calibration and validation
#' @param y the response vector.
#' @param y0 the reponse vector without noise \code{sigmay}.
#' @param sigmay the uncertainty on \code{y}.
#' @param sigmaondes the standard deviation of the Gaussians.
#' @author François Wahl Louna Alsouki
#' @seealso `browseVignettes("dual.spls")`
#'
#' @examples
#' ### load dual.spls library
#' library(dual.spls)
#' ### parameters
#' n <- 100
#' p <- 50
#' nondes <- 20
#' sigmaondes <- 0.5
#' data.benchmark=BCHMK(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
#'
#' X <- data.benchmark$X
#' y <- data.benchmark$y
#'
#' ###plotting the data
#' plot(X[1,],type='l',ylim=c(0,max(X)),main='Benchmark data', ylab='X',col=1)
#' for (i in 2:n){
#'  lines(X[i,],col=i)
#'}
#' @export


#############################
# This function divides the y in cells of equal range
# and attributes a type to the observations according to the corresponding cell
# then divides the data into calibration validation using the function FW_LAS_f1

# Arguments:
# X : dataset to split in calibration and validation
# pcal : proportion to put in calibration
# Datatype : number of the cell to which each observation belongs
# y : values of the responses
# ncells : number of cells to build if necessary

# Values:
# indcal : index of the observations to put in the calibration dataset
# indval : index of the observations to put in the validation dataset
#############################

FW_LAS_calval<- function(X,pcal,Datatype=NULL,y=NULL,ncells=10)
{
  if (is.null(Datatype) & is.null(y)){
    stop('Error in FW_LAS_calval: if Datatype=NULL, y should not be NULL' )
  }

  # type of the observed values
  if (is.null(Datatype)){
    Datatype=FW_LAS_type(y,ncells)
  }

  # nb elts in calibration for each cell
  ycounts=sapply(1:ncells,function(u) sum(Datatype==u) )
  Listecal=ceiling(ycounts*pcal/100)

  # index of calibration/validation
  indcal=FW_LAS_f1(X,Datatype,Listecal)
  indval=1:dim(X)[1]
  indval=indval[-indcal]

  return(list(indcal=indcal, indval=indval))
}
