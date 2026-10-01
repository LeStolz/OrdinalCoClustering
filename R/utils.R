#' Calculate the complete log-likelihood of a model.
#'
#' @param model A `model` object.
#' @param ... Additional arguments.
#'
#' @return The complete log-likelihood.
#' @export
log_likelihood <- function(model, ...) {
  UseMethod("log_likelihood")
}


#' Calculate the ICL of a model.
#'
#' @param model A `model` object.
#' @param ... Additional arguments.
#'
#' @return The ICL value.
#' @export
icl <- function(model, ...) {
  UseMethod("icl")
}