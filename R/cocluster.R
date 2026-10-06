#' Co-cluster data by fitting a latent block model using any observation model
#' and SEM-Gibbs estimation.
#'
#' This function performs simultaneous clustering of the rows and columns of an
#' ordinal data matrix. It relies on the Latent Block Model (LBM) and uses a
#' Stochastic Expectation-Maximization (SEM) algorithm coupled with Gibbs sampling.
#'
#' @param x Matrix of ordinal observations. Missing values (`NA`) are natively supported.
#' @param K Initial number of row clusters.
#' @param L Initial number of column clusters.
#' @param model Observation model. This should be a list of functions defining the 
#'   probability distribution of a block. For ordinal data, use `cub(m)` where `m` 
#'   is the maximum number of ordinal categories.
#' @param n_init Number of random initializations. Multiple initializations are 
#'   strongly recommended to avoid local maxima.
#' @param n_iterations Number of iterations for the SEM-Gibbs algorithm.
#' @param n_gibbs_iterations Number of Gibbs iterations for the SEM-Gibbs algorithm 
#'   (typically performed within the E-step).
#' @param eps Convergence threshold for the SEM-Gibbs algorithm.
#' @param burn_in Number of burn-in iterations for the SEM-Gibbs algorithm. 
#'   Parameter estimates are averaged after this burn-in period.
#' @param choose_n_clusters Logical. If `TRUE`, the function will greedily search 
#'   the neighborhood of `(K, L)` and select the optimal number of clusters based 
#'   on the Integrated Completed Likelihood (ICL) criterion.
#'
#' @return An object of class `cocluster` containing:
#' \describe{
#'   \item{lbms}{List of all models evaluated.}
#'   \item{selected_lbm}{The optimal model selected based on ICL.}
#'   \item{selected_lbm$icl}{ICL value of the selected model.}
#'   \item{selected_lbm$p$row}{Estimated global probabilities of row-cluster membership.}
#'   \item{selected_lbm$p$col}{Estimated global probabilities of column-cluster membership.}
#'   \item{selected_lbm$z$row}{Final row-cluster assignments.}
#'   \item{selected_lbm$z$col}{Final column-cluster assignments.}
#'   \item{selected_lbm$block_parameters}{Estimated CUB parameters (`pi` and `xi`) of each block.}
#' }
#' 
#' @examples
#' # 1. Create a small random ordinal matrix (20 rows, 10 columns, 1 to 5 scale)
#' set.seed(42)
#' x <- matrix(sample(1:5, 200, replace = TRUE), nrow = 20, ncol = 10)
#' 
#' # 2. Fit the co-clustering model 
#' # (using small iterations for the sake of the example)
#' result <- cocluster(
#'   x = x, 
#'   K = 2, 
#'   L = 2, 
#'   model = cub(5), 
#'   n_init = 2, 
#'   n_iterations = 10, 
#'   n_gibbs_iterations = 2,
#'   choose_n_clusters = FALSE
#' )
#' 
#' # 3. Inspect the optimal row and column assignments
#' result$selected_lbm$z$row
#' result$selected_lbm$z$col
#'
#' @export
cocluster <- function(x,
                      K,
                      L,
                      model = cub,
                      n_init = 20,
                      n_iterations = 100,
                      n_gibbs_iterations = 10,
                      eps = 1e-9,
                      burn_in = floor(n_iterations / 2),
                      choose_n_clusters = TRUE) {

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (K < 1 || K > nrow(x)) {
    stop("K must be between 1 and the number of rows.")
  }

  if (L < 1 || L > ncol(x)) {
    stop("L must be between 1 and the number of columns.")
  }

  if (!is.list(model)) {
    stop("model must be a list.")
  }

  required_functions <- c(
    "probability",
    "log_probability",
    "estimate_parameters",
	  "number_of_parameters"
  )

  missing_functions <- required_functions[
    !required_functions %in% names(model)
  ]

  if (length(missing_functions) > 0) {
    stop(
      "Model is missing: ",
      paste(missing_functions, collapse = ", ")
    )
  }

  fit_result <- fit_models(
    x = x,
    K = K,
    L = L,
    model = model,
    n_init = n_init,
    n_iterations = n_iterations,
    n_gibbs_iterations = n_gibbs_iterations,
    eps = eps,
    burn_in = burn_in
  )
  lbms <- c(fit_result$lbms)
  selected_lbm <- fit_result$selected_lbm

  # https://inria.hal.science/hal-01658589/document
  #
  if (choose_n_clusters) {
    repeat {
      prev_selected_lbm <- selected_lbm
      neighbors <- get_neighbors(prev_selected_lbm$K$row, prev_selected_lbm$K$col, nrow(x), ncol(x))

      for (neighbor in neighbors) {
        fit_result <- fit_models(
          x = x,
          K = neighbor[1],
          L = neighbor[2],
          model = model,
          n_init = n_init,
          n_iterations = n_iterations,
          n_gibbs_iterations = n_gibbs_iterations,
          eps = eps,
          burn_in = burn_in
        )

        lbms <- c(lbms, fit_result$lbms)

        if (is.null(selected_lbm) || fit_result$selected_lbm$icl > selected_lbm$icl) {
          selected_lbm <- fit_result$selected_lbm
        }
      }

      if (identical(selected_lbm, prev_selected_lbm)) {
        break
      }
    }
  }

  structure(
    list(
      lbms = lbms,
      selected_lbm = selected_lbm
    ),
    class = "cocluster"
  )
}


fit_models <- function(x,
                       K,
                       L,
                       model,
                       n_init,
                       n_iterations,
                       n_gibbs_iterations,
                       eps,
                       burn_in) {
  lbms <- vector("list", n_init)
  selected_lbm <- NULL

  for (init in seq_len(n_init)) {
    fit <- lbm(
      x = x,
      K = K,
      L = L,
      model = model,
      n_iterations = n_iterations,
      n_gibbs_iterations = n_gibbs_iterations,
      eps = eps,
      burn_in = burn_in
    )

    if (is.null(selected_lbm) || fit$icl > selected_lbm$icl) {
      selected_lbm <- fit
    }
	  lbms[[init]] <- fit
  }

  list(
    lbms = lbms,
    selected_lbm = selected_lbm
  )
}


get_neighbors <- function(K, L, N, D) {
  neighbors <- list()

  if (K < N) {
    neighbors[[length(neighbors) + 1]] <- c(K + 1, L)
  }

  if (K > 1) {
    neighbors[[length(neighbors) + 1]] <- c(K - 1, L)
  }

  if (L < D) {
    neighbors[[length(neighbors) + 1]] <- c(K, L + 1)
  }

  if (L > 1) {
    neighbors[[length(neighbors) + 1]] <- c(K, L - 1)
  }

  neighbors
}