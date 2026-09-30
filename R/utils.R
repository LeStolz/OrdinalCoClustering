#' Calculate the complete log-likelihood
#'
#' @param x A `model` object.
#' @param ... Additional arguments.
#'
#' @return The complete log-likelihood.
#' @export
log_likelihood <- function(x, ...) {
  UseMethod("log_likelihood")
}


#' Calculate the ICL
#'
#' @param x A `model` object.
#' @param ... Additional arguments.
#'
#' @return The ICL value.
#' @export
icl <- function(x, ...) {
  UseMethod("icl")
}