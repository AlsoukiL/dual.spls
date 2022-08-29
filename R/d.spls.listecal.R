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
  l=rep(0,ncells)
  i=0
  rep=floor(n*pcal/100)
  while (i < rep)
  {
    for (j in 1:ncells)
    {
      if (ycounts[j]>0 && i<rep)
      {
        ycounts[j]=ycounts[j]-1
        l[j]=l[j]+1
        i=i+1
      }
    }

  }

  Listecal=l
  return(Listecal)
}

