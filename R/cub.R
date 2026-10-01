#' Probability of a response under the CUB model
#'
#' Mixture of a deliberate choice (shifted binomial, parameter `xi`) and
#' a random choice (uniform). Note: `xi` close to 0 means high ratings.
#'
#' @param r Vector of responses between 1 and m.
#' @param m Number of ordinal categories.
#' @param pi_ Probability of a deliberate choice, in (0, 1).
#' @param xi Feeling parameter, in (0, 1).
#' @return A vector of probabilities, same length as `r`.
#' @noRd
cub_probability <- function(r, m, pi_, xi) {
  # P(R = r) = pi * b_r(xi) + (1 - pi) / m
  # dbinom(r - 1, m - 1, 1 - xi) = probabilité d'avoir r - 1 succès sur
  # m - 1 essais, avec une probabilité de succès 1 - xi : c'est la
  # binomiale décalée (partie "choix réfléchi").
  # (1 - pi_) / m : partie "hasard", loi uniforme sur 1..m.
  pi_ * dbinom(r - 1, size = m - 1, prob = 1 - xi) + (1 - pi_) / m
}

#' Log-probability of each cell of a block under the CUB model
#'
#' @param block Matrix (or vector) of responses between 1 and m.
#' @param m Number of ordinal categories.
#' @param pi_ Probability of a deliberate choice.
#' @param xi Feeling parameter.
#' @return Same shape as `block`, with the log-probability of each response.
#'   Missing responses (NA) give 0 so that they are ignored in a sum.
#' @noRd
cub_log_probability <- function(block, m, pi_, xi) {
  # On travaille en log : un produit de nombreuses probabilités tombe à 0
  # sur ordinateur, alors qu'une somme de logs reste stable.
  log_p <- log(cub_probability(block, m, pi_, xi))

  # Une réponse manquante n'apporte pas d'information : on met 0 pour
  # que sum() ne renvoie pas NA.
  log_p[is.na(block)] <- 0
  log_p
}

#' E-step of the EM algorithm for the CUB model
#'
#' Posterior probability that each response comes from the deliberate
#' choice component rather than from the uniform component.
#'
#' @param r Vector of responses between 1 and m (no NA).
#' @param m Number of ordinal categories.
#' @param pi_ Current value of pi.
#' @param xi Current value of xi.
#' @return A vector tau with values in 0, 1, same length as `r`.
#' @noRd
cub_e_step <- function(r, m, pi_, xi) {
  # Règle de Bayes : tau_i = pi * b(r_i) / P(r_i)
  # = part du "choix réfléchi" dans la probabilité totale de la note r_i.
  # tau_i est l'espérance de la variable latente z_i (1 = choix réfléchi).
  pi_ * dbinom(r - 1, size = m - 1, prob = 1 - xi) /
    cub_probability(r, m, pi_, xi)
}

#' M-step of the EM algorithm for the CUB model
#'
#' @param r Vector of responses between 1 and m (no NA).
#' @param tau Posterior probabilities from the E-step.
#' @param m Number of ordinal categories.
#' @return A list with the updated `pi_` and `xi`.
#' @noRd
cub_m_step <- function(r, tau, m) {
  # pi : proportion moyenne de réponses issues du choix réfléchi.
  pi_new <- mean(tau)

  # xi : maximise sum(tau * log b(r)). On dérive et on annule, ce qui donne
  # 1 - xi = moyenne pondérée de (r - 1) divisée par (m - 1).
  # Si sum(tau) = 0 (division par 0), xi n'est pas identifiable : valeur neutre.
  if (sum(tau) == 0) {
    xi_new <- 0.5
  } else {
    xi_new <- 1 - sum(tau * (r - 1)) / ((m - 1) * sum(tau))
  }

  list(pi_ = pi_new, xi = xi_new)
}

#' EM algorithm for the CUB model on a vector of responses
#'
#' @param r Vector of responses between 1 and m (may contain NA).
#' @param m Number of ordinal categories.
#' @param max_iterations Maximum number of EM iterations.
#' @param tolerance Convergence threshold on the change of the parameters.
#' @param bound Keeps pi and xi inside bound, 1 - bound to avoid log(0).
#' @return A named list `list(pi_ = , xi = )`.
#' @noRd
cub_em <- function(r, m, max_iterations = 200, tolerance = 1e-8,
                   bound = 1e-3) {
  r <- r[!is.na(r)]

  # Bloc vide (cluster vide pendant le SEM) : paramètres neutres,
  # sinon mean(tau) donnerait NaN.
  if (length(r) == 0) {
    return(list(pi_ = 0.5, xi = 0.5))
  }

  # Initialisation neutre. L'EM converge vers un maximum local :
  # le résultat peut dépendre de ce point de départ.
  pi_ <- 0.5
  xi <- 0.5

  for (iteration in seq_len(max_iterations)) {
    tau <- cub_e_step(r, m, pi_, xi)          # étape E
    updated <- cub_m_step(r, tau, m)          # étape M

    # On borne pi et xi dans [bound, 1 - bound] : évite log(0) = -Inf
    # dans la log-vraisemblance (blocs presque constants).
    new_pi <- min(max(updated$pi_, bound), 1 - bound)
    new_xi <- min(max(updated$xi, bound), 1 - bound)

    # Critère d'arrêt : les paramètres ne bougent presque plus.
    change <- abs(new_pi - pi_) + abs(new_xi - xi)
    pi_ <- new_pi
    xi <- new_xi

    if (change < tolerance) break
  }

  list(pi_ = pi_, xi = xi)
}

#' Create a CUB observation model for the latent block model
#'
#' Returns the list of functions expected by [cocluster()]:
#' `probability`, `log_probability`, `estimate_parameters`
#' and `number_of_parameters`.
#'
#' @param m Number of ordinal categories (responses are in 1, ..., m).
#' @return A list defining the CUB model.
#' @export
make_cub_model <- function(m) {
  # Fonction "usine" : m est mémorisé par les 3 fonctions ci-dessous
  # (fermeture), car un bloc peut ne pas contenir toutes les notes.
  list(
    probability = function(x, parameters) {
      cub_probability(x, m, parameters$pi_, parameters$xi)
    },
    log_probability = function(block, parameters) {
      cub_log_probability(block, m, parameters$pi_, parameters$xi)
    },
    # Appelée par le SEM à chaque étape M, une fois par bloc.
    estimate_parameters = function(block) {
      cub_em(as.vector(block), m)
    },
    # pi et xi : utilisé dans la pénalité de l'ICL.
    number_of_parameters = 2
  )
}