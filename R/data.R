#' Example dataset for TDVS package demonstration
#'
#' A simulated dataset containing response and predictor variables for testing
#' thresholded variable selection.
#'
#' @format A named list with 2 components:
#' \describe{
#'   \item{y}{Numeric response matrix (100 rows, 1 column)}
#'   \item{X}{Numeric predictor matrix (100 rows, 8 columns)}
#' }
#' @usage data(data_tdvs)
#' @examples
#' data(data_tdvs)
#' str(data_tdvs)
#' plot(data_tdvs$X[,1], data_tdvs$y)
#' @name data_tdvs
#' @docType data
NULL
