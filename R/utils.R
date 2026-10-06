#' Calculate the complete log-likelihood of a model.
#'
#' @param model An object of class `model` containing:
#' \describe{
#'   \item{x}{The data matrix.}
#'   \item{z}{A list containing the cluster assignments.}
#'   \item{p}{A list containing the prior probabilities of cluster membership.}
#'   \item{block_parameters}{Estimated observation model parameters of each cluster.}
#'   \item{model}{The observation model containing the ``log_probability`` function.}
#' }
#' @param ... Additional arguments.
#'
#' @return The complete log-likelihood.
#' @export
log_likelihood <- function(model, ...) {
  UseMethod("log_likelihood")
}


#' Calculate the ICL of a model.
#'
#' @param model An object of class `model` containing:
#' \describe{
#'   \item{N}{The number of observations.}
#'   \item{K}{The number of clusters.}
#'   \item{x}{The data matrix.}
#'   \item{z}{A list containing the cluster assignments.}
#'   \item{p}{A list containing the prior probabilities of cluster membership.}
#'   \item{block_parameters}{Estimated observation model parameters of each cluster.}
#'   \item{model}{The observation model containing the ``number_of_parameters`` and the ``log_probability`` function.}
#' }
#' @param ... Additional arguments.
#'
#' @return The ICL value.
#' @export
icl <- function(model, ...) {
  UseMethod("icl")
}