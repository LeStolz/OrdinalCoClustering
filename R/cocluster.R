#' Co-cluster data by fitting a Latent Block Model using a provided observation model,
#' SEM-Gibbs estimation and model selection using ICL.
#'
#' @param x Matrix of observations.
#' @param K Initial number of row clusters. If a vector is provided, a grid search is performed.
#' @param L Initial number of column clusters. If a vector is provided, a grid search is performed.
#' @param model Observation model as a list of functions defining the probability distribution of a block containing:
#' \describe{
#'   \item{probability}{Function to compute the probability of given responses under the obseration model.}
#'   \item{log_probability}{Function to compute the log probability of given responses under the obseration model.}
#'   \item{estimate_parameters}{Function to estimate the parameters of the observation model.}
#'   \item{number_of_parameters}{The number of parameters in the observation model.}
#' }
#' For ordinal data, use `cub(m)` where `m` is the number of ordinal categories.
#' @param n_init Number of random initializations to avoid local maxima.
#' @param n_iterations Number of iterations for the SEM-Gibbs algorithm.
#' @param n_gibbs_iterations Number of Gibbs iterations for the SEM-Gibbs algorithm in the E step.
#' @param eps Convergence threshold for the SEM-Gibbs algorithm.
#' @param burn_in Number of burn-in iterations for the SEM-Gibbs algorithm.
#'   Parameter estimates are averaged after this burn-in period.
#' @param choose_n_clusters If `TRUE`, the function will greedily search (https://inria.hal.science/hal-01658589/document)
#'   the neighborhood of `(K, L)` and select the optimal number of clusters based
#'   on the Integrated Completed Likelihood (ICL) criterion.
#'
#' @return An object of class `cocluster` containing:
#' \describe{
#'   \item{lbms}{List of all \code{\link{lbm}} models evaluated.}
#'   \item{selected_lbm}{An object of class \code{\link{lbm}} corresponding to the optimal model selected based on ICL.}
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
#'   n_init = 1,
#'   n_iterations = 2,
#'   n_gibbs_iterations = 1,
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
                      model,
                      n_init = 20,
                      n_iterations = 100,
                      n_gibbs_iterations = 10,
                      eps = 1e-9,
                      burn_in = floor(n_iterations / 2),
                      choose_n_clusters = TRUE) {

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (any(K < 1) || any(K > nrow(x))) {
    stop("K must be between 1 and the number of rows.")
  }

  if (any(L < 1) || any(L > ncol(x))) {
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

  if (length(K) > 1 || length(L) > 1) {
    if (choose_n_clusters) {
      warning("choose_n_clusters is ignored when multiple values are provided for K or L (grid search is performed instead).")
    }

    grid <- expand.grid(K = K, L = L)
    lbms <- list()
    selected_lbm <- NULL

    for (i in seq_len(nrow(grid))) {
      fit_result <- fit_models(
        x = x,
        K = grid$K[i],
        L = grid$L[i],
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
  } else {
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
  }

  structure(
    list(
      lbms = lbms,
      selected_lbm = selected_lbm
    ),
    class = "cocluster"
  )
}


#' Fit multiple Latent Block Models and select the best one
#'
#' @param x Matrix of observations.
#' @param K Number of row clusters.
#' @param L Number of column clusters.
#' @param model Observation model.
#' @param n_init Number of random initializations.
#' @param n_iterations Number of iterations for the SEM-Gibbs algorithm.
#' @param n_gibbs_iterations Number of Gibbs iterations for the SEM-Gibbs algorithm.
#' @param eps Convergence threshold for the SEM-Gibbs algorithm.
#' @param burn_in Number of burn-in iterations for the SEM-Gibbs algorithm.
#'
#' @return A list containing:
#' \describe{
#'   \item{lbms}{List of all \code{\link{lbm}} models evaluated.}
#'   \item{selected_lbm}{The optimal \code{\link{lbm}} model selected based on ICL.}
#' }
#' @noRd
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


#' Get neighboring dimensions for ICL greedy search
#'
#' @param K Number of row clusters.
#' @param L Number of column clusters.
#' @param N Number of rows in the data matrix.
#' @param D Number of columns in the data matrix.
#'
#' @return A list of vectors `c(K_new, L_new)` representing the valid neighboring cluster dimensions.
#' @noRd
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