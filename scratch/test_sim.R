# CUB simulation function
simulate_cub <- function(n, m, pi_, xi) {
  deliberate <- rbinom(n, size = m - 1, prob = 1 - xi) + 1
  random <- sample(1:m, size = n, replace = TRUE)
  is_deliberate <- rbinom(n, size = 1, prob = pi_) == 1
  ifelse(is_deliberate, deliberate, random)
}

set.seed(42)
N <- 100
D <- 80
K_true <- 3
L_true <- 2
m <- 5

p_row <- c(0.4, 0.3, 0.3)
p_col <- c(0.6, 0.4)

z_row <- sample(1:K_true, size = N, replace = TRUE, prob = p_row)
z_col <- sample(1:L_true, size = D, replace = TRUE, prob = p_col)

pi_true <- matrix(c(0.9, 0.8, 0.2, 0.9, 0.8, 0.8), nrow = K_true, ncol = L_true, byrow = TRUE)
xi_true <- matrix(c(0.1, 0.9, 0.5, 0.5, 0.8, 0.2), nrow = K_true, ncol = L_true, byrow = TRUE)

x <- matrix(NA, nrow = N, ncol = D)
for (i in 1:N) {
  for (j in 1:D) {
    k <- z_row[i]
    l <- z_col[j]
    x[i, j] <- simulate_cub(1, m, pi_true[k, l], xi_true[k, l])
  }
}

source('R/cub.R')
source('R/lbm.R')
source('R/utils.R')
source('R/cocluster.R')

result <- cocluster(x = x, K = 3, L = 2, model = cub(m), n_init = 2, n_iterations = 20, n_gibbs_iterations = 2, choose_n_clusters = FALSE)

print(result$selected_lbm$p$row)
