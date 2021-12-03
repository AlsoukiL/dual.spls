#' Splits data into calibration and validation sets using the FWLAS splitting method
#' @description
#' The function \code{d.spls.FWLAS.calval} divides the data \code{X} into a calibration and a validation. According to the values
#' of the response vector, the observations are divided into groups. The Method uses
#' the Kennard and Stone startegy for each group at a time and according to the number of calibration desired
#' from each group, it selects the calibration points.
#' @usage d.spls.FWLAS.calval(X,pcal=NULL,Datatype=NULL,y=NULL,ncells=10,Listecal=NULL)
#' @param X a numeric matrix of predictors values.
#' @param pcal a positive integer between 0 and 100. \code{pcal} is the percentage
#' of calibration samples to be selected. Default value is NULL, meaning as long as \code{Listecal} is
#' specified, \code{pcal} is not necessary.
#' @param Datatype A vecor of index specifying each observation belonging to wich group index.
#' Default value is \code{NULL}, meaning the function will use the internal function \code{type} to compute the vector for \code{ncells}.
#' If \code{NULL}, parameter \code{y} should be specified.
#' @param y a numeric vector of responses. Default value is \code{NULL}, meaning as long as \code{Datatype} is specified,
#' \code{y} is not necessary.
#' @param ncells a positive integer. \code{ncells} is the number of groups dividing the observations. If
#' \code{Datatype} is not specified, the function devides the observations into \code{ncells} groups. Default value
#' is \code{10}.
#' @param Listecal a numeric vector specifying how many observations from each group should be selected as calibration.
#' Default value is \code{NULL}, meaning the function will consider a percentage of \code{pcal} from each group
#' to be in the calibration set. If \code{NULL}, parameter \code{pcal} should be specified.
#' @details
#' The modified Kennard and Stone algorithm allows to select samples using the classical Kennard and stone on
#' each group of observations one by one. It starts by selecting the point that is the furthest away from the centroid.
#' This point is assigned as the calibration set and is removed from the list of candidates. Then, it identifies to which
#' group belongs this first observation and consider the group \eqn{g} that comes after.
#' It computes the distance \eqn{\delta_{P_{i,g}}} between the remaining points
#' \eqn{P_{i,g}} belonging to the group the group \eqn{g} and the calibration point assigned. The point with the
#' largest \eqn{\delta_{P_{i,g}}} is selected, removed from the set and the procedure moves on to the group that comes
#' after.
#'
#' When there is more than one calibration sample, the procedure computes the distance between each \eqn{P_{i,g}} from
#' the concerned group and each \eqn{P_{i,cal}} from the calibration set. The selected candidate verify the following distance:
#' \deqn{d_{selected}=\max\limits_{P_{i,g}} ( \min\limits_{P_{i,cal}} (d_{P_{i,g},P_{i,cal}}) )}
#'
#' Once each of the vector \code{Listecal} elements are null, the procedure is done.
#' @return A \code{list} of the following attributes
#' \item{indcal}{a numeric vector giving the row indices of the input data selected for calibration.}
#' \item{indval}{a numeric vector giving the row indices of the remaining observations.}
#' @author Louna Alsouki François Wahl
#' @seealso [dual.spls::d.spls.calval()],[dual.spls::d.spls.type()],`browseVignettes("dual.spls")`
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
#' ###calibration parameters for split1
#' pcal <- 70
#' ncells <- 3
#'
#' split1 <- d.spls.FWLAS.calval(X=X,pcal=pcal,y=y,ncells=ncells)
#'
#' ###plotting split1
#' plot(X[split1$indcal,1],X[split1$indcal,2],xlab="Variable 1",
#' ylab="Variable 2",pch=19,col="red",main="Calibration and validation split1")
#' points(X[split1$indval,1],X[split1$indval,2],pch=19,col="green")
#' legend("topright", legend = c("Calibration points", "Validation points"),
#' cex = 0.8, col = c("red","green"), pch = c(19,19))
#'
#' ###calibration parameters for split2
#' ncells <- 3
#' dimtype=floor(n/3)
#' Datatype <- c(rep(1,dimtype),rep(2,dimtype),rep(3,(n-dimtype*2)))
#'
#' L1=floor(0.7*length(which(Datatype==1)))
#' L2=floor(0.8*length(which(Datatype==2)))
#' L3=floor(0.6*length(which(Datatype==3)))
#' Listecal <- c(L1,L2,L3)
#'
#' split2 <- d.spls.FWLAS.calval(X=X,y=y,Datatype=Datatype,Listecal=Listecal)
#'
#' ###plotting split2
#' plot(X[split2$indcal,1],X[split2$indcal,2],xlab="Variable 1",
#' ylab="Variable 2",pch=19,col="red",main="Calibration and validation split2")
#' points(X[split2$indval,1],X[split2$indval,2],pch=19,col="green")
#' legend("topright", legend = c("Calibration points", "Validation points"),
#'  cex = 0.8, col = c("red","green"), pch = c(19,19))
#'
#' @export


d.spls.FWLAS.calval<- function(X,pcal=NULL,Datatype=NULL,y=NULL,ncells=10,Listecal=NULL)
{
  if (is.null(Datatype) & is.null(y)){
    stop('if Datatype=NULL, y should not be NULL' )
  }

  # type of the observed values
  if (is.null(Datatype)){
    Datatype=d.spls.type(y,ncells)
  }

  if (is.null(Listecal) & is.null(pcal)){
    stop('if Listecal=NULL, pcal should not be NULL' )
  }

  #percentage of calibration
  if(is.null(Listecal)){
    # nb elts in calibration for each cell
    ycounts=sapply(1:ncells,function(u) sum(Datatype==u) )
    Listecal=ceiling(ycounts*pcal/100)
  }

  if(max(Datatype) != length(Listecal)){
    stop('length of Listecal does not match with values of Datatype' )
  }

  # index of calibration/validation
  indcal=d.spls.calval(X=X,Xtype=Datatype,Listecal=Listecal)
  indval=1:dim(X)[1]
  indval=indval[-indcal]

  return(list(indcal=indcal, indval=indval))
}
