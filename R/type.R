#' Splitting the observations into groups
#' @keywords internal
#' @description
#' The internal function \code{type} divides the response vector \code{y} in \code{cells}
#' of equal range and attributes a index type to the observations according to the corresponding cell
#' @param y a vector.
#' @param ncells a positive integer. \code{ncells} is the number of subsamples desired.
#' #' @return A vecor of index specifying each observation belonging to wich group index.
#' @author François Wahl Louna Alsouki
#' @seealso `browseVignettes("dual.spls")`


type<- function(y,ncells)
{
  # ybreaks
  ymM=max(y)-min(y)
  ybreaks=seq(min(y)-1/(2*ncells)*ymM,max(y)+1/(2*ncells)*ymM,length.out=ncells+1) #

  # Affecting each class to X
  Datatype=rep(0,length(y))
  for (ib in 1:ncells){
    sel= (y>=ybreaks[ib]) & (y<ybreaks[ib+1])
    Datatype[sel]=ib
  }
  Datatype[which.max(y)]=ib
  Datatype=as.vector(Datatype)

  return(Datatype)
}
