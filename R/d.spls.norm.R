#' Computing the l1 and l2 norm
#' @keywords internal
#' @description
#' The internal functions \code{norm1} and \code{norm2} computes the l1 and l2 norm of any vector
#' @usage d.spls.norm1(x)
#' @usage d.spls.norm2(x)
#' @param x a numeric vector.
#' @return A positive integer of the norm value.
#' @author Louna Alsouki François Wahl
#' @seealso `browseVignettes("dual.spls")`


d.spls.norm1 <- function(x) {
  d <- sum(abs(x))
  return(d)
}


d.spls.norm2 <- function(x) {
  d <- sqrt(sum(x^2))
  return(d)
}
