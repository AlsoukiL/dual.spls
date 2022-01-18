#' Splits data into calibration and validation sets according to wich group belongs each observation
#' @keywords internal
#' @description
#' The function \code{d.spls.calval} divides the data \code{X} into a calibration and a validation set using
#' the Kennard and Stone strategy for each group at a time and according to the number of calibration desired
#' from each group.
#' @param X a numeric matrix.
#' @param Xtype a vector of index specifying to which group belongs each observation.
#' @param Listecal a vector specifying how many observations from each group should be selected as calibration.
#' @return a numeric vector giving the row indices of the input data selected for calibration
#' @author Louna Alsouki François Wahl
#' @seealso [dual.spls::d.spls.type()],[dual.spls::d.spls.FWLAS.calval()],`browseVignettes("dual.spls")`
#' @importFrom pdist pdist
#'
d.spls.calval<- function(X,Xtype,Listecal)
{

  if (length(Xtype)< dim(X)[1] || all.equal(sort(unique(Xtype)),c(1:max(Xtype))) != T)
  {
    stop("some observations are badly indexed")
  }
  for (i in 1:max(Xtype))
  {
  if (length(which(Xtype==i)) < Listecal[i] )
  {
    stop("number of calibration observations chosen exceeds the overall number for the group ", i)
  }
}
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

  # centroid of X
  G=apply(X,2,mean)

  # distance matrix between each observation of X and the centroid
  D_XG=as.matrix(pdist::pdist(X,Y=G))
  # the furthest observation from the centroid
  N=which.max(D_XG)

  # first element of indcal
  ical=1
  indcal[ical]=N

  # the group to which it belongs
  Listecal[Xtype[N]]=Listecal[Xtype[N]]-1

  # the counter index
  i=Xtype[indcal[1]]
  while (sum(Listecal)>0)
  {
    # finding which group to consider next
    for (i1 in ((i+1):(i+n_grp))){
      if (i1 <= n_grp)
        igr=i1
      else
        igr=i1-n_grp

      if (Listecal[igr]>0)
        break
    }
    i=igr

    # find the maxmin point for calibration
    indreste=1:n_exp
    indreste=indreste[-indcal]
    indreste=indreste[Xtype[indreste]==igr]
    D_calreste=as.matrix(pdist::pdist(X[indcal,],X[indreste,]))
    D_calrestemin=apply(D_calreste,2,min)
    indnew=indreste[which.max(D_calrestemin)]

    ical=ical+1

    indcal[ical]=indnew
    Listecal[igr]=Listecal[igr]-1

  }
  return(indcal)
}
