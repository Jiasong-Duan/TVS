#' Example dataset for TVS package demonstration
#'
#' A simulated dataset containing response and predictor variables for testing thresholded variable selection.
#'
#' @format A named list with 2 components:
#' \describe{
#'   \item{y}{Numeric response matrix (100 rows, 1 column)}
#'   \item{X}{Numeric predictor matrix (100 rows, 8 columns)}
#' }
#' @usage data(data_tvs)
#' @examples
#' data(data_tvs)
#' str(data_tvs)
#' plot(data_tvs$X[,1], data_tvs$y)
"data_tvs"
