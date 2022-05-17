#' Determine the number of observations to select from each cell
#' @keywords internal
#' @description
#' The function \code{d.spls.listecal} detemines the number of observations to select as calibration from
#' each cell.
#' @param Xtype a vector of index specifying to which group belongs each observation.
#' @param pcal a positive integer between 0 and 100. \code{pcal} is the percentage
#' of calibration samples to be selected.
#' @return a numeric vector specifying how many observations from each group should be selected as calibration.
#' @author Louna Alsouki François Wahl
#' @seealso [dual.spls::d.spls.type],[dual.spls::d.spls.calval],[dual.spls::d.spls.split]
#' @importFrom pdist pdist

d.spls.listecal<- function(Xtype,pcal)
{
  ycounts=table(Xtype)
  ncells=length(ycounts)
  n=length(Xtype)
  igr=which(ycounts < n/ncells)
  not.igr=1:ncells
  if (sum(igr) > 0) not.igr=not.igr[-igr]
  Listecal=rep(0,ncells)
  Listecal[igr]=ycounts[igr]
  n.igr=sum(Listecal[igr])
  pcal2=100*(pcal/100*n-n.igr)/(n-n.igr)
  Listecal[not.igr]=ceiling(ycounts[not.igr]*pcal2/100)

  return(Listecal)
}
