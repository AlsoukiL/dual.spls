#' Splits data into calibration and validation sets according to wich group belongs each observation
#' @keywords internal
#' @description
#' The function \code{modified.KS} divides the data \code{X} into a calibration and a validation set using
#' the Kennard and Stone startegy for each group at a time and according to the number of calibration desired
#' from each group
#' @param X a numeric matrix.
#' @param Xtype a vector of index specifying to wich group belongs each observation.
#' @param Listecal a vector specifying how many observations from each group should be selected as calibration.
#' @return a numeric vector giving the row indices of the input data selected for calibration
#' @seealso [dual.spls::type()],[dual.spls::calval()],`browseVignettes("dual.spls")`
#' @importFrom pdist pdist
#'
modified.KS<- function(X,Xtype,Listecal)
{

  nref=sum(Listecal)

  ###################################
  # Sorting the data
  ###################################
  # groups
  grp=unique(Xtype)
  # number of groups
  n_grp=length(grp)
  # number of observations
  n_exp=dim(X)[1]

  # number of calibration points desired
  ncal=sum(Listecal)


  # intializing the vector of calibration index
  indcal=array(0,sum(Listecal))

  ###################################

  # Centroid of X
  G=apply(X,2,mean)

  # Distance matrix between each observation of X and the centroid
  D_XG=as.matrix(pdist(X,G))
  # the furthest observation from the centroid
  N=which.max(D_XG)

  # first element of indcal
  ical=1
  indcal[ical]=N
  # The group to which it belongs
  Listecal[Xtype[N]]=Listecal[Xtype[N]]-1

  #the counter index
  i=Xtype[indcal[1]]
  while (sum(Listecal)>0)
  {
    # Finding which group to consider next
    for (i1 in ((i+1):(i+n_grp))){
      if (i1 <= n_grp)
        igr=i1
      else
        igr=i1-n_grp

      if (Listecal[igr]>0)
        break
    }
    i=igr

    # Find the maxmin point for calibration
    indreste=1:n_exp
    indreste=indreste[-indcal]
    indreste=indreste[Xtype[indreste]==igr]
    D_calreste=as.matrix(pdist(X[indcal,],X[indreste,]))
    D_calrestemin=apply(D_calreste,2,min)
    indnew=indreste[which.max(D_calrestemin)]

    ical=ical+1

    indcal[ical]=indnew
    Listecal[igr]=Listecal[igr]-1

  }
  return(indcal)
}
