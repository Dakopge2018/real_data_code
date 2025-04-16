# Libraries
library(depmixS4)
library(ggplot2)
library(MASS) # Pour générer des échantillons bivariés normaux

# Model parameters
set.seed(123)
n <- 365 * 100
states <- 2
degree_obs_pol <- 1
degree_trans_pol <- 1
period <- 365
time <- 1:n
n_sim <- 1000

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
  trig_covs <- data.frame(time = time, trend = floor(time/period))
  for (d in 1:degree) {
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * time / period)
    trig_covs[[paste0("sin_", d)]] <- sin(2 * pi * d * time / period)
  }
  return(trig_covs)
}

# Simulate two Normal observations
# Simulate two independent Normal observations
simulate_observations <- function(states, n, period, mu, sigma, delta_normal, degree, trend_params_normal) {
  obs_normal1 <- numeric(n)
  obs_normal2 <- numeric(n)
  
  for (t in 1:n) {
    k <- states[t]
    
    # Mean for the first normal variable
    mean_t1 <- mu[k, 1] + trend_params_normal[k, 1] * floor(t / period) + 
      sum(delta_normal[k, , 1] * c(sin(2 * pi * (1:degree) * t / period), cos(2 * pi * (1:degree) * t / period)))
    
    # Mean for the second normal variable
    mean_t2 <- mu[k, 2] + trend_params_normal[k, 2] * floor(t / period) + 
      sum(delta_normal[k, , 2] * c(sin(2 * pi * (1:degree) * t / period), cos(2 * pi * (1:degree) * t / period)))
    
    # Variances for the two variables (independent, so no covariance)
    var1 <- sigma[[k]][1, 1]
    var2 <- sigma[[k]][2, 2]
    
    # Generate independent normal samples
    obs_normal1[t] <- rnorm(1, mean = mean_t1, sd = sqrt(var1))
    obs_normal2[t] <- rnorm(1, mean = mean_t2, sd = sqrt(var2))
  }
  
  return(list(normal1 = obs_normal1, normal2 = obs_normal2))
}

# Parameters
beta <- list(
  c(1, 0.7, 0.5),
  c(-1, -0.6, 0.7)
)
beta_matrix <- matrix(beta, nrow = states, ncol = states - 1, byrow = TRUE)

# Emission parameters
# Emission parameters
mu <- matrix(c(0.2, 2, 0.04, 3), nrow = 2, byrow = TRUE)  # Moyennes pour chaque état
sigma <- list(
  matrix(c(1, 0, 0, 3), nrow = 2),  # Variances pour l'état 1
  matrix(c(1.5, 0, 0, 2), nrow = 2)  # Variances pour l'état 2
)
delta_normal <- array(c(0.5, 0.3, -0.4, 0.2, 0.1, -0.2, 0.3, 0.4), dim = c(2, 2, 2))
trend_params_normal <- matrix(c(0.005, 0.01, 0.002, 0.008), nrow = 2, byrow = TRUE)

# Simulate data
simulated_states <- simulate_states(n, period, beta_matrix, states, degree_trans_pol)
observations <- simulate_observations(simulated_states, n, period, mu, sigma, delta_normal, degree_obs_pol, trend_params_normal)
trig_covs <- generate_trig_covariates(time, period, degree_trans_pol)
# Access the simulated variables
head(observations$normal1)  # Première variable normale
head(observations$normal2)  # Deuxième variable normale

data <- data.frame(
  obs_normal1 = observations$normal1,
  obs_normal2 = observations$normal2,
  trig_covs
)

# Fit model with two Normal responses
transition_formula <- as.formula(paste("~", paste(names(trig_covs)[c(-1,-2)], collapse = " + ")))
mod <- depmix(
  list(obs_normal1 ~ sin_1 + cos_1 + trend, obs_normal2 ~ sin_1 + cos_1 + trend),
  data = data,
  nstates = 2,
  family = list(gaussian(), gaussian()),
  transition = transition_formula,
)

set.seed(1)
fitted_model <- multistart(mod,
  nstarts = 10,  # 10 initialisations
  initIters = 10,  # 10 itérations EM pour chaque initialisation
  emcontrol = em.control(
    maxit = 500,  # Max 500 itérations EM
    tol = 1e-08,  # Tolérance pour convergence
    crit = "relative",  # Critère de convergence
    random.start = TRUE  # Randomisation des initialisations
  )
)

summary(fitted_model, which = 'transition')
summary(fitted_model, which = 'response')

# Graphiques pour comparer les observations simulées et ajustées
fitted_states <- posterior(fitted_model, type = 'viterbi')$state

# comparison_normal1 <- data.frame(
#   Time = 1:n,
#   States = simulated_states,
#   Simulated = observations$normal1,
#   Fitted = fitted_model@response[[1]][[1]]@parameters$fitted
# )

# comparison_normal2 <- data.frame(
#   Time = 1:n,
#   States = simulated_states,
#   Simulated = observations$normal2,
#   Fitted = fitted_model@response[[2]][[1]]@parameters$fitted
# )

# # Graphique pour la première normale
# ggplot(comparison_normal1[1:1000,], aes(x = Time)) +
#   geom_line(aes(y = Simulated, color = "Simulated"), size = 1) +
#   geom_line(aes(y = Fitted, color = "Fitted"), size = 1, linetype = "dashed") +
#   labs(title = "Comparison of Simulated and Fitted Observations (Normal 1)",
#        subtitle = "First 1000 Time Points",
#        y = "Value",
#        color = "Legend") +
#   theme_minimal()

# # Graphique pour la deuxième normale
# ggplot(comparison_normal2[1:1000,], aes(x = Time)) +
#   geom_line(aes(y = Simulated, color = "Simulated"), size = 1) +
#   geom_line(aes(y = Fitted, color = "Fitted"), size = 1, linetype = "dashed") +
#   labs(title = "Comparison of Simulated and Fitted Observations (Normal 2)",
#        subtitle = "First 1000 Time Points",
#        y = "Value",
#        color = "Legend") +
#   theme_minimal()
# Fonction pour calculer les probabilités de transition simulées
compute_simulated_transition_probs <- function(t, period, beta_matrix, K, degree) {
  P <- matrix(0, nrow = K, ncol = K)
  
  for (i in 1:K) {
  
    for (j in 2:K) {  # Éviter la diagonale principale
      P[i, j] <- beta_matrix[[i, j - 1]][[1]][[1]][[1]]  # Constante
      for (d in 1:degree) {  # Ajouter les termes trigonométriques
        P[i, j] <- P[i, j] + 
          beta_matrix[[i, j - 1]][[1]][[2 * d]][[1]] * sin(2 * pi * d * t / period) + 
          beta_matrix[[i, j - 1]][[1]][[2 * d +1]][[1]] * cos(2 * pi * d * t / period)
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
  Time = 1:n_sim,
  Simulated = simulated_Q11[1:n_sim],
  Estimated = estimated_Q11[1:n_sim]
)

comparison_Q22 <- data.frame(
  Time = 1:n_sim,
  Simulated = simulated_Q22[1:n_sim],
  Estimated = estimated_Q22[1:n_sim]
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


# Modified simulation from fitted model for Poisson
simulate_observations_from_fitted <- function(fitted_model, states, n, period, degree) {
  obs_poisson <- numeric(n)
  obs_normal <- numeric(n)
  for (t in 1:n) {
    # if (states[t] == 1) {
    #   k <- 2
    # } else {
    #   k <- 1
    # }
    k <- states[t]
    response_model_poisson <- fitted_model@response[[k]][[1]]@parameters
    response_model_normal <- fitted_model@response[[k]][[2]]@parameters
    response_params_poisson <- response_model_poisson$coefficients
    response_params_normal <- response_model_normal$coefficients
       # Poisson
    lambda_t <- exp(response_params_poisson[1] + response_params_poisson[4] * floor(t/period)) * exp(sum(response_params_poisson[2]* sin(2 * pi * (1:degree) * t/period), response_params_poisson[3]*cos(2 * pi * (1:degree) * t/period))) # nolint: line_length_linter.
    obs_poisson[t] <- rpois(1, lambda = lambda_t)
  
    # Normal
    mean_t <- response_params_normal[1] + response_params_normal[4] * floor(t/period) + sum(response_params_normal[2]* sin(2 * pi * (1:degree) * t /period), response_params_normal[3]*cos(2 * pi * (1:degree) * t/period))
    obs_normal[t] <- rnorm(1, mean = mean_t, sd = response_model_normal$sd)
  }
  return(list(poisson = obs_poisson, normal = obs_normal))
}

# Rest of comparison and plotting code remains the same
# fitted_states <- simulate_states_from_fitted(fitted_model, n, period, degree_trans_pol)
fitted_states = posterior(fitted_model, type ='viterbi')$state
fitted_observations <- simulate_observations_from_fitted(fitted_model, fitted_states, n, period, degree_obs_pol)


comparison_poisson <- data.frame(
 Time = 1:n,
 States = simulated_states,
 Simulated = observations$poisson,
 Fitted = fitted_observations$poisson
)

comparison_normal <- data.frame(
 Time = 1:n,
 States = simulated_states,
 Simulated = observations$normal,
 Fitted = fitted_observations$normal
)
mse <- mean((comparison_normal$Simulated - comparison_normal$Fitted)^2)
mae <- mean(abs(comparison_normal$Simulated - comparison_normal$Fitted))
cor_coef <- cor(comparison_normal$Simulated, comparison_normal$Fitted)

mse <- mean((comparison_poisson$Simulated - comparison_poisson$Fitted)^2)
mae <- mean(abs(comparison_poisson$Simulated - comparison_poisson$Fitted))
cor_coef <- cor(comparison_poisson$Simulated, comparison_poisson$Fitted)

# Graphique pour les observations Poisson
ggplot(comparison_poisson[1:100,], aes(x = Time)) +
  geom_line(aes(y = Simulated, color = "Simulated"), size = 1) +
  geom_line(aes(y = Fitted, color = "Fitted"), size = 1, linetype = "dashed") +
  labs(title = "Comparison of Simulated and Fitted Observations (Poisson)",
       subtitle = "First 1000 Time Points",
       y = "Count",
       color = "Legend") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_color_manual(values = c("Simulated" = "blue", "Fitted" = "red"))

# Graphique pour les observations normales
ggplot(comparison_normal[1:1000,], aes(x = Time)) +
  geom_line(aes(y = Simulated, color = "Simulated"), size = 1) +
  geom_line(aes(y = Fitted, color = "Fitted"), size = 1, linetype = "dashed") +
  labs(title = "Comparison of Simulated and Fitted Observations (Normal)",
       subtitle = "First 1000 Time Points",
       y = "Value",
       color = "Legend") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_color_manual(values = c("Simulated" = "blue", "Fitted" = "red")