#' Create a Latent Block Model
#'
#' @param x Matrix of observations.
#' @param K Number of row clusters.
#' @param L Number of column clusters.
#' @param model Observation model.
#'
#' @return An object of class `lbm`.
#' @export
lbm <- function(x,
                K,
                L,
                model,
                n_iterations = 100,
                n_gibbs_iterations = 10,
                eps = 1e-3,
                burn_in = floor(n_iterations / 2)) {

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
    burn_in = floor(n_iterations / 2)
  )
}


#' Initialize a Latent Block Model
#'
#' @param x Matrix of observations.
#' @param K Number of row clusters.
#' @param L Number of column clusters.
#' @param model Observation model.
#'
#' @return An object of class `lbm`.
initialize_lbm <- function(x,
                           K,
                           L,
                           model) {
  N <- c(row = nrow(x), col = ncol(x))
  K <- c(row = K, col = L)

  # https://www.sciencedirect.com/science/article/pii/S0167947325000647
  # Random initialization making sure every cluster has at least one observation
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


#' Estimate parameters for each LBM block
#'
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


#' Calculate the complete log-likelihood
#'
#' @param lbm An `lbm` object.
#'
#' @return The complete log-likelihood.
#' @export
log_likelihood.lbm <- function(lbm) {
  x <- lbm$x
  z <- lbm$z
  block_parameters <- lbm$block_parameters
  model <- lbm$model

  log_likelihood <- sum(log(lbm$p$row[z$row])) + sum(log(lbm$p$col[z$col]))

  for (k in seq_len(lbm$K$row)) {
    for (l in seq_len(lbm$K$col)) {
      row_indices <- which(z$row == k)
      col_indices <- which(z$col == l)

      if (length(row_indices) == 0 || length(col_indices) == 0) {
        next
      }

      block <- x[row_indices, col_indices, drop = FALSE]

      block_log_p <- model$log_probability(block, block_parameters[[k, l]])

      log_likelihood <- log_likelihood + sum(block_log_p)
    }
  }

  log_likelihood
}


#' Calculate the ICL
#'
#' @param lbm An `lbm` object.
#'
#' @return The ICL value.
#' @export
icl.lbm <- function(lbm) {
  N <- lbm$N$row
  D <- lbm$N$col
  K <- lbm$K$row
  L <- lbm$K$col
  v <- lbm$model$number_of_parameters

  log_likelihood_value <- log_likelihood(lbm)

  penalty <- (
    (K - 1) * log(N) +
    (L - 1) * log(D) +
    K * L * v * log(N * D)
  ) / 2

  log_likelihood_value - penalty
}


#' Calculate conditional cluster probabilities
#'
#' @noRd
calculate_p_zi_ks__x_theta <- function(lbm, dimension, index) {
  K <- lbm$K[dimension]
  p_zi_ks__x_theta <- numeric(K)

  for (k in seq_len(K)) {
    p_zi_k__x_theta <- lbm$p[[dimension]][k]

    if (dimension == "row") {
      for (j in seq_len(lbm$N$col)) {
        l <- lbm$z$col[j]
        x <- lbm$x[index, j]
        parameters <- lbm$block_parameters[[k, l]]

        p_zi_k__x_theta <- p_zi_k__x_theta * lbm$model$probability(x, parameters)
      }
    } else if (dimension == "col") {
      for (i in seq_len(lbm$N$row)) {
        l <- lbm$z$row[i]
        x <- lbm$x[i, index]
        parameters <- lbm$block_parameters[[k, l]]

        p_zi_k__x_theta <- p_zi_k__x_theta * lbm$model$probability(x, parameters)
      }
    } else {
      stop("dimension must be either 'row' or 'col'.")
    }

    p_zi_ks__x_theta[k] <- p_zi_k__x_theta
  }

  p_zi_ks__x_theta / sum(p_zi_ks__x_theta)
}


#' Sampling cluster assignments using Gibbs
#'
#' @noRd
gibbs_sample_cluster <- function(lbm, dimension, hard_assignment = FALSE) {
  N <- lbm$N[dimension]
  K <- lbm$K[dimension]

  for (index in seq_len(N)) {
    p_zi_ks__x_theta <- calculate_p_zi_ks__x_theta(
      lbm = lbm,
      dimension = dimension,
      index = index,
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


#' SEM-Gibbs algorithm
#' https://inria.hal.science/inria-00494796/document
#'
#' @param lbm An `lbm` object.
#' @param n_gibbs_iterations Number of Gibbs iterations in the E-S step.
#' @param n_iterations Number of SEM iterations.
#' @param eps Convergence threshold.
#' @param burn_in Number of iterations to discard as burn-in.
#'
#' @noRd
sem_gibbs <- function(lbm, n_iterations = 100, n_gibbs_iterations = 10, eps = 1e-3, burn_in = floor(n_iterations / 2)) {
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
  }
  else {
    parameter_history <- parameter_history[seq.int(n_iter_completed, n_iter_completed)]
  }
  lbm <- calculate_final_parameters(lbm, parameter_history, n_iterations)

  lbm
}


#' Average parameter estimates obtained after burn-in
#'
#' @param lbm An `lbm` object.
#' @param parameter_history List of parameter estimates.
#'
#' @return Averaged model parameters.
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
            parameter_history, function(x) x$block_parameters[[k, l]][[name]], numeric(1)
          ))
        }
	    )

      names(block_parameters[[k, l]]) <- names(first_block_parameters[[k, l]])
    }
  }

  lbm$p <- p
  lbm$block_parameters <- block_parameters

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