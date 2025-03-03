# Libraries
library(depmixS4)
library(ggplot2)

# Model parameters
set.seed(123)
n <- 365 * 5
states <- 2
degree_obs_pol <- 1
degree_trans_pol <- 1
period <- 365
time <- 1:n

# Function to compute the seasonal transition matrix
compute_transition_matrix <- function(t, period, beta_matrix, K, degree) {
  P <- matrix(0, nrow = K, ncol = K)
  
  for (i in 1:K) {
    for (j in 2:K) {  # Avoid diagonal
      P[i, j] <- beta_matrix[[i, j - 1]][1]  # Constant term
      for (d in 1:degree) {  # Add seasonal components
        P[i, j] <- P[i, j] + 
          beta_matrix[[i, j - 1]][2 * d] * sin(2 * pi * d * t / period) + 
          beta_matrix[[i, j - 1]][2 * d + 1] * cos(2 * pi * d * t / period)
      }
    }
  }
  
  # Convert P to probability matrix Q
  Q <- matrix(0, nrow = K, ncol = K)
  for (i in 1:K) {
    exp_values <- exp(P[i, 2:K])
    Q[i, 2:K] <- exp_values / (1 + sum(exp_values))
    Q[i, 1] <- 1 / (1 + sum(exp_values))
  }
  return(Q)
}

# Function to simulate states
simulate_states <- function(n, period, beta_matrix, K, degree) {
  states <- numeric(n)
  states[1] <- sample(1:K, 1)
  
  for (t in 2:n) {
    Q_t <- compute_transition_matrix(t, period, beta_matrix, K, degree)
    states[t] <- sample(1:K, size = 1, prob = Q_t[states[t - 1], ])
  }
  return(states)
}

# Generate trigonometric covariates
generate_trig_covariates <- function(time, period, degree) {
  trig_covs <- data.frame(time = time, trend = time)
  for (d in 1:degree) {
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * time / period)
    trig_covs[[paste0("sin_", d)]] <- sin(2 * pi * d * time / period)
  }
  return(trig_covs)
}

# Simulate Poisson and Normal observations
simulate_observations <- function(states, n, period, mu, nu, delta_poisson, delta_normal, sd, degree, trend_params) {
  obs_poisson <- numeric(n)
  obs_normal <- numeric(n)
  
  for (t in 1:n) {
    k <- states[t]
    
    # Poisson
    lambda_t <- exp(mu[k] + trend_params[k] * t) * exp(sum(delta_poisson[k, ] * c(sin(2 * pi * (1:degree) * t / period), cos(2 * pi * (1:degree) * t / period))))
    obs_poisson[t] <- rpois(1, lambda = lambda_t)
    
    # Normal
    mean_t <- nu[k] + trend_params[k] * t + sum(delta_normal[k, ] * c(sin(2 * pi * (1:degree) * t / period), cos(2 * pi * (1:degree) * t / period)))
    obs_normal[t] <- rnorm(1, mean = mean_t, sd = sd[k])
  }
  return(list(poisson = obs_poisson, normal = obs_normal))
}

# Parameters
beta <- list(
  c(1, 0.7, 0.5),
  c(-1, -0.6, 0.7)
)
beta_matrix <- matrix(beta, nrow = states, ncol = states - 1, byrow = TRUE)

# Emission parameters
mu <- c(2, 0.7)
nu <- c(0, 2)
delta_poisson <- matrix(c(1, 0.7, -0.9, 0.1), nrow = states, byrow = TRUE)
delta_normal <- matrix(c(0.5, 0.3, -0.4, 0.2), nrow = states, byrow = TRUE)
sd <- c(1, 0.5)
trend_params <- c(0.001, -0.001)  # Specific trend for each state

# Simulate data
simulated_states <- simulate_states(n, period, beta_matrix, states, degree_trans_pol)
trig_covs <- generate_trig_covariates(time, period, degree_trans_pol)
observations <- simulate_observations(simulated_states, n, period, mu, nu, delta_poisson, delta_normal, sd, degree_obs_pol, trend_params)

data <- data.frame(
  obs_poisson = observations$poisson,
  obs_normal = observations$normal,
  trig_covs
)

# Fit model with Poisson response
transition_formula <- as.formula(paste("~", paste(names(trig_covs)[c(-1,-2)], collapse = " + ")))
mod <- depmix(
  list(obs_poisson ~ sin_1 + cos_1 + trend, obs_normal ~ sin_1 + cos_1 + trend),
  data = data,
  nstates = 2,
  family = list(poisson(), gaussian()),
  transition = transition_formula,
)

set.seed(1)
fitted_model  <- multistart(mod,
  nstarts = 10,  # 10 initialisations
  initIters = 10,  # 10 itérations EM pour chaque initialisation
  emcontrol = em.control(
    maxit = 500,  # Max 500 itérations EM
    tol = 1e-08,  # Tolérance pour convergence
    crit = "relative",  # Critère de convergence
    random.start = TRUE,  # Randomisation des initialisations
    classification = "Hard"  # Classification 
    ))
summary(fitted_model, which = 'transition')
summary(fitted_model, which = 'response')
# Fonction pour calculer les probabilités de transition simulées
compute_simulated_transition_probs <- function(t, period, beta_matrix, K, degree) {
  P <- matrix(0, nrow = K, ncol = K)
  
  for (i in 1:K) {
    for (j in 2:K) {  # Éviter la diagonale principale
      P[i, j] <- beta_matrix[[i, j - 1]][[1]][[1]][[1]]  # Constante
      for (d in 1:degree) {  # Ajouter les termes trigonométriques
        P[i, j] <- P[i, j] + 
          beta_matrix[[i, j - 1]][[1]][[2 * d]][[1]] * cos(2 * pi * d * t / period) + 
          beta_matrix[[i, j - 1]][[1]][[2 * d +1]][[1]] * sin(2 * pi * d * t / period)
      }
    }
  }
  
  # Transformer P en une matrice de probabilités Q
  Q <- matrix(0, nrow = K, ncol = K)
  for (i in 1:K) {
    exp_values <- exp(P[i, 2:K])  # Exponentielles des valeurs
    Q[i, 2:K] <- exp_values / (1 + sum(exp_values))
    Q[i, 1] <- 1 / (1 + sum(exp_values))  # Proba pour rester
  }
  return(Q)
}

# Calculer les probabilités de transition simulées
simulated_Q11 <- sapply(1:n, function(t) compute_transition_matrix(t, period, beta_matrix, states, degree_trans_pol)[1, 1])
simulated_Q22 <- sapply(1:n, function(t) compute_transition_matrix(t, period, beta_matrix, states, degree_trans_pol)[2, 2])

# Extraire les probabilités de transition estimées
estimated_beta <- list(
    list(fitted_model@transition[[1]]@parameters$coefficients[,2]),
    list(fitted_model@transition[[2]]@parameters$coefficients[,2])
    )

est_beta_matrix <- matrix(estimated_beta, nrow = states, ncol = states - 1, byrow = TRUE)
# Comparer les probabilités de transition simulées et estimées
estimated_Q11 <- sapply(1:n, function(t) compute_simulated_transition_probs(t, period, est_beta_matrix, states, degree_trans_pol)[1, 1])
estimated_Q22 <- sapply(1:n, function(t) compute_simulated_transition_probs(t, period, est_beta_matrix, states, degree_trans_pol)[2, 2])
comparison_Q11 <- data.frame(
  Time = 1:200,
  Simulated = simulated_Q11[1:200],
  Estimated = estimated_Q11[1:200]
)

comparison_Q22 <- data.frame(
  Time = 1:200,
  Simulated = simulated_Q22[1:200],
  Estimated = estimated_Q22[1:200]
)

# Tracer les probabilités de transition pour comparaison
ggplot(comparison_Q11, aes(x = Time)) +
  geom_line(aes(y = Simulated, color = "Simulated")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "Comparison of Simulated and Estimated Q11(t)",
       y = "Q11(t)",
       color = "Legend") +
  theme_minimal()
  ggplot(comparison_Q22, aes(x = Time)) +
  geom_line(aes(y = Simulated, color = "Simulated")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "Comparison of Simulated and Estimated Q22(t)",
       y = "Q22(t)",
       color = "Legend") +
  theme_minimal()
