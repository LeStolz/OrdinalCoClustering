#' Create a Latent Block Model to estimate and model a matrix of observations.
#'
#' @param x Matrix of observations.
#' @param K Number of row clusters.
#' @param L Number of column clusters.
#' @param model Observation model as a list of functions defining the probability distribution of a block containing:
#' \describe{
#'   \item{probability}{Function to compute the probability of given responses under the observation model.}
#'   \item{log_probability}{Function to compute the log probability of given responses under the observation model.}
#'   \item{estimate_parameters}{Function to estimate the parameters of the observation model.}
#'   \item{number_of_parameters}{The number of parameters in the observation model.}
#' }
#' @param n_iterations Number of iterations for the SEM-Gibbs algorithm used to estimate the model.
#' @param n_gibbs_iterations Number of Gibbs iterations for the SEM-Gibbs algorithm used to estimate the model.
#' @param eps Convergence threshold for the SEM-Gibbs algorithm.
#' @param burn_in Number of burn-in iterations for the SEM-Gibbs algorithm.
#'
#' @return An object of class `lbm` containing:
#' \describe{
#'   \item{x}{The data matrix.}
#'   \item{N}{A list containing the number of rows (`row`) and columns (`col`) of the data matrix.}
#'   \item{K}{A list containing the number of row clusters (`row`) and column clusters (`col`).}
#'   \item{z}{A list containing the row-cluster (`row`) and column-cluster (`col`) assignments.}
#'   \item{p}{A list containing the prior probabilities of row-cluster (`row`) and column-cluster (`col`) membership.}
#'   \item{block_parameters}{A matrix of the estimated observation model parameters for each block, indexed by `[[row_cluster, col_cluster]]`.}
#'   \item{model}{The observation model used.}
#'   \item{icl}{The Integrated Completed Likelihood (ICL) value of the model.}
#' }
#' @export
lbm <- function(x,
                K,
                L,
                model,
                n_iterations,
                n_gibbs_iterations,
                eps,
                burn_in) {

  initial_lbm <- initialize_lbm(
    x = x,
    K = K,
    L = L,
    model = model
  )

  sem_gibbs(
    lbm = initial_lbm,
    n_iterations = n_iterations,
    n_gibbs_iterations = n_gibbs_iterations,
    eps = eps,
    burn_in = burn_in
  )
}


#' Initialize a Latent Block Model
#'
#' @param x Matrix of observations.
#' @param K Number of row clusters.
#' @param L Number of column clusters.
#' @param model Observation model as a list of functions containing the `estimate_parameters` function.
#'
#' @return An object of class \code{\link{lbm}} containing the initial random assignments and parameters.
#' @noRd
initialize_lbm <- function(x,
                           K,
                           L,
                           model) {
  N <- list(row = nrow(x), col = ncol(x))
  K <- list(row = K, col = L)

  # Initialisation aléatoire en s'assurant que chaque cluster est non vide
  repeat {
    z <- list(
      row = sample(seq_len(K$row), size = N$row, replace = TRUE),
      col = sample(seq_len(K$col), size = N$col, replace = TRUE)
    )

    p <- list(
      row = tabulate(z$row, nbins = K$row) / N$row,
      col = tabulate(z$col, nbins = K$col) / N$col
    )

    if (all(p$row > 0) && all(p$col > 0)) {
      break
    }
  }

  lbm <- structure(
    list(
      x = x,
      N = N,
      K = K,
      z = z,
      p = p,
      block_parameters = NULL,
      model = model,
      icl = NULL
    ),
    class = "lbm"
  )

  lbm$block_parameters <- estimate_block_parameters(lbm = lbm)

  lbm
}


#' Estimate parameters for each Latent Block Model block
#'
#' @param lbm An object of class \code{\link{lbm}}.
#'
#' @return A matrix containing the estimated parameters for each block.
#' @noRd
estimate_block_parameters <- function(lbm) {
  x <- lbm$x
  z <- lbm$z
  K <- lbm$K
  model <- lbm$model

  parameters <- matrix(
    vector("list", K$row * K$col),
    nrow = K$row,
    ncol = K$col
  )

  for (k in seq_len(K$row)) {
    for (l in seq_len(K$col)) {
      row_indices <- which(z$row == k)
      col_indices <- which(z$col == l)

      block <- x[row_indices, col_indices, drop = FALSE]

      parameters[[k, l]] <- model$estimate_parameters(block)
    }
  }

  parameters
}


#' Calculate the complete log-likelihood of the Latent Block Model
#'
#' @param model An object of class \code{\link{lbm}}.
#'
#' @return The complete log-likelihood.
#' @export
log_likelihood.lbm <- function(model) {
  x <- model$x
  z <- model$z
  block_parameters <- model$block_parameters
  obs_model <- model$model

  log_likelihood <- sum(log(model$p$row[z$row])) + sum(log(model$p$col[z$col]))

  for (k in seq_len(model$K$row)) {
    for (l in seq_len(model$K$col)) {
      row_indices <- which(z$row == k)
      col_indices <- which(z$col == l)

      if (length(row_indices) == 0 || length(col_indices) == 0) {
        next
      }

      block <- x[row_indices, col_indices, drop = FALSE]

      block_log_p <- obs_model$log_probability(block, block_parameters[[k, l]])

      log_likelihood <- log_likelihood + sum(block_log_p)
    }
  }

  log_likelihood
}


#' Calculate the ICL of the Latent Block Model
#'
#' @param model An object of class \code{\link{lbm}}.
#'
#' @return The ICL value.
#' @export
icl.lbm <- function(model) {
  N <- model$N$row
  D <- model$N$col
  K <- model$K$row
  L <- model$K$col
  v <- model$model$number_of_parameters

  log_likelihood_value <- log_likelihood(model)

  penalty <- (
    (K - 1) * log(N) +
    (L - 1) * log(D) +
    K * L * v * log(N * D)
  ) / 2

  log_likelihood_value - penalty
}


#' Calculate the conditional cluster probabilities
#' that a row or column belongs to each cluster given the other dimension's assignments and parameters.
#'
#' The probabilities are calculated in log-space to avoid numerical underflow.
#'
#' @param lbm An object of class \code{\link{lbm}}.
#' @param dimension String of either `"row"` or `"col"` to specify if the index refers to a row or column
#' @param index The index of the row or column to calculate probabilities for.
#'
#' @return A vector of conditional probabilities corresponding to each cluster.
#' @noRd
calculate_p_zi_ks__x_theta <- function(lbm, dimension, index) {
  K <- lbm$K[[dimension]]
  log_p_zi_ks <- numeric(K)

  for (k in seq_len(K)) {
    if (lbm$p[[dimension]][k] == 0) {
      log_p_zi_ks[k] <- -Inf
      next
    }

    log_p <- log(lbm$p[[dimension]][k])

    if (dimension == "row") {
      for (l in seq_len(lbm$K$col)) {
        col_indices <- which(lbm$z$col == l)
        if (length(col_indices) == 0) next

        block <- lbm$x[index, col_indices, drop = FALSE]
        parameters <- lbm$block_parameters[[k, l]]

        log_p <- log_p + sum(lbm$model$log_probability(block, parameters))
      }
    } else if (dimension == "col") {
      for (l in seq_len(lbm$K$row)) {
        row_indices <- which(lbm$z$row == l)
        if (length(row_indices) == 0) next

        block <- lbm$x[row_indices, index, drop = FALSE]
        parameters <- lbm$block_parameters[[l, k]]

        log_p <- log_p + sum(lbm$model$log_probability(block, parameters))
      }
    } else {
      stop("dimension must be either 'row' or 'col'.")
    }

    log_p_zi_ks[k] <- log_p
  }

  max_log_p <- max(log_p_zi_ks)
  if (is.infinite(max_log_p)) {
    return(rep(1 / K, K))
  }

  p_zi_ks <- exp(log_p_zi_ks - max_log_p)
  p_zi_ks / sum(p_zi_ks)
}


#' Sample cluster assignments using Gibbs
#'
#' Gibbs calculates the conditional probabilities of cluster assignments for each row or column,
#' and samples the new cluster assignments from these probabilities.
#'
#' @param lbm An object of class \code{\link{lbm}}.
#' @param dimension String of either `"row"` or `"col"`, specifying which dimension to sample.
#' @param hard_assignment If `TRUE`, assigns to the cluster with maximum probability instead of sampling.
#'
#' @return The updated \code{\link{lbm}} object with new cluster assignments.
#' @noRd
gibbs_sample_cluster <- function(lbm, dimension, hard_assignment = FALSE) {
  N <- lbm$N[[dimension]]
  K <- lbm$K[[dimension]]

  for (index in seq_len(N)) {
    p_zi_ks__x_theta <- calculate_p_zi_ks__x_theta(
      lbm = lbm,
      dimension = dimension,
      index = index
    )

    if (hard_assignment) {
      new_zi <- which.max(p_zi_ks__x_theta)
    } else {
      new_zi <- sample(
        seq_len(K),
        size = 1,
        prob = p_zi_ks__x_theta
      )
    }

    lbm$z[[dimension]][index] <- new_zi
  }

  lbm
}


#' Perform the SEM-Gibbs algorithm
#' https://inria.hal.science/inria-00494796/document
#'
#' The SEM-Gibbs algorithm iteratively samples cluster assignments
#' and estimates parameters for the Latent Block Model using maximum likelihood estimation.
#'
#' @param lbm An object of class \code{\link{lbm}}.
#' @param n_iterations Number of SEM iterations.
#' @param n_gibbs_iterations Number of Gibbs iterations in the E step.
#' @param eps Convergence threshold.
#' @param burn_in Number of iterations to discard as burn-in.
#'
#' @return An object of class \code{\link{lbm}} containing the final estimated parameters and assignments.
#' @noRd
sem_gibbs <- function(lbm, n_iterations, n_gibbs_iterations, eps, burn_in) {
  parameter_history <- vector("list", n_iterations)
  prev_log_likelihood <- -Inf
  n_iter_completed <- 0

  for (iteration in seq_len(n_iterations)) {
    # SE
    for (gibbs_iteration in seq_len(n_gibbs_iterations)) {
      lbm <- gibbs_sample_cluster(
        lbm = lbm,
        dimension = "row"
      )

      lbm <- gibbs_sample_cluster(
        lbm = lbm,
        dimension = "col"
      )
    }

    # M
    lbm$p <- list(
      row = tabulate(lbm$z$row, nbins = lbm$K$row) / lbm$N$row,
      col = tabulate(lbm$z$col, nbins = lbm$K$col) / lbm$N$col
    )

    lbm$block_parameters <- estimate_block_parameters(lbm = lbm)

    parameter_history[[iteration]] <- list(
      p = lbm$p,
      block_parameters = lbm$block_parameters
    )
    n_iter_completed <- n_iter_completed + 1

    cur_log_likelihood <- log_likelihood(lbm)
    if (abs(cur_log_likelihood - prev_log_likelihood) < eps) {
      break
    }
    prev_log_likelihood <- cur_log_likelihood
  }

  if (burn_in < n_iter_completed) {
    parameter_history <- parameter_history[seq.int(burn_in + 1, n_iter_completed)]
  } else {
    parameter_history <- parameter_history[seq.int(n_iter_completed, n_iter_completed)]
  }

  lbm <- calculate_final_parameters(lbm, parameter_history, n_iterations)
}


#' Average parameter estimates obtained after burn-in
#'
#' After burn-in, The SEM-Gibbs algorithm samples parameters around the maximum likelihood estimates.
#' This function averages the estimates after burn-in to provide a stable estimate of the model parameters
#' and calculates the final cluster assignments based on these averaged parameters.
#'
#' @param lbm An object of class \code{\link{lbm}}.
#' @param parameter_history List of parameter estimates.
#' @param n_iterations Number of iterations for the SEM-Gibbs algorithm.
#'
#' @return An object of class \code{\link{lbm}} containing the averaged parameters and final cluster assignments.
#'
#' @noRd
calculate_final_parameters <- function(lbm, parameter_history, n_iterations) {
  first_block_parameters <- parameter_history[[1]]$block_parameters
  block_parameters <- first_block_parameters

  N <- length(parameter_history)

  p <- list(
    row = Reduce("+", lapply(parameter_history, function(x) x$p$row)) / N,
    col = Reduce("+", lapply(parameter_history, function(x) x$p$col)) / N
  )

  K <- nrow(first_block_parameters)
  L <- ncol(first_block_parameters)

  for (k in seq_len(K)) {
    for (l in seq_len(L)) {
      block_parameters[[k, l]] <- lapply(
        names(first_block_parameters[[k, l]]),
        function(name) {
          mean(vapply(
            parameter_history,
            function(x) x$block_parameters[[k, l]][[name]],
            numeric(1)
          ))
        }
      )

      names(block_parameters[[k, l]]) <- names(first_block_parameters[[k, l]])
    }
  }

  lbm$p <- p
  lbm$block_parameters <- block_parameters

  # Calculate the final cluster assignments based on the averaged parameters after burn-in
  for (iteration in seq_len(n_iterations)) {
    old_z <- lbm$z

    lbm <- gibbs_sample_cluster(
      lbm = lbm,
      dimension = "row",
      hard_assignment = TRUE
    )
    lbm <- gibbs_sample_cluster(
      lbm = lbm,
      dimension = "col",
      hard_assignment = TRUE
    )

    if (identical(old_z, lbm$z)) {
      break
    }
  }

  lbm$icl <- icl(lbm)

  lbm
}