#' Dual Sparse Partial Least Squares (Dual-SPLS) Near Infrared data
#' @description
#' This dataset contains absorbance spectra of refined petroleum samples measured between 4500 nm and 9000 nm interval. The
#' oil refinery is an industrial process plant where petroleum is transformed and refined into useful products such as gasoline,
#' diesel fuel, asphalt base, fuel oils, heating oil, kerosene... Each product is distinguishable by its characteristics like
#' its density also available in the dataset.
#' @usage data(d.spls.NIR)
#' @format A \code{data.frame} of 208 observations. The NIR data contains 1557 variables and the density is a one dimensional vector.
#' @details
#' The NIR data is composed of 1557 variables that represent a regular sequence of the interval \eqn{[4500,9000]}.
#'
#' This data is useful when building a Dual-SPLS regression with \code{NIR} data as the predictor and \code{density} as a response.
"d.spls.NIR"
