#' Determine if object is scalar.
#'
#' Simple checker if the value is 1-dimensional
#' @param x a vector of length 1
#'
#' @return A logical TRUE or FALSE
#' @export
#'
#' @examples
#' is.scalar(1)
#' is.scalar(c(1,2))
#' @importFrom plyr ddply
is.scalar <- function(x) {
  # plyr::ddply()
  is.atomic(x) && length(x) == 1L
}
