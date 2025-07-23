# Libraries
library(depmixS4)
library(ggplot2)

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

# Simulate Poisson and Normal observations
simulate_observations <- function(states, n, period, mu, nu, delta_poisson, delta_normal, sd, degree, trend_params_poisson, trend_params_normal) {
  obs_poisson <- numeric(n)
  obs_normal <- numeric(n)
  
  for (t in 1:n) {
    k <- states[t]
    
    # Poisson
    lambda_t <- exp(mu[k] + trend_params_poisson[k] * floor(t/period)) * exp(sum(delta_poisson[k, ] * c(sin(2 * pi * (1:degree) * t / period), cos(2 * pi * (1:degree) * t / period))))
    obs_poisson[t] <- rpois(1, lambda = lambda_t)
    
    # Normal
    mean_t <- nu[k] + trend_params_normal[k] * floor(t/period) + sum(delta_normal[k, ] * c(sin(2 * pi * (1:degree) * t / period), cos(2 * pi * (1:degree) * t / period)))
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
mu <- c(0.1, 0.05)
nu <- c(0, 2)
delta_poisson <- matrix(c(1, 0.7, -0.9, 0.1), nrow = states, byrow = TRUE)
delta_normal <- matrix(c(0.5, 0.3, -0.4, 0.2), nrow = states, byrow = TRUE)
sd <- c(1, 0.5)
trend_params_poisson <- c(0.001, 0.0012)# Specific trend for each state
trend_params_normal <- c(0.005, 0.01) # Specific trend for each state
# Simulate data
simulated_states <- simulate_states(n, period, beta_matrix, states, degree_trans_pol)
trig_covs <- generate_trig_covariates(time, period, degree_trans_pol)
observations <- simulate_observations(simulated_states, n, period, mu, nu, delta_poisson, delta_normal, sd, degree_obs_pol, trend_params_poisson, trend_params_normal)

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
  family = list(poisson(link = "log"), gaussian()),
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
    random.start = TRUE  # Randomisation des initialisations
 #   classification = "Hard"  # Classification 
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
  Time = 1:n_sim,
  Simulated = simulated_Q11[1:n_sim],
  Estimated = estimated_Q22[1:n_sim]
)

comparison_Q22 <- data.frame(
  Time = 1:n_sim,
  Simulated = simulated_Q22[1:n_sim],
  Estimated = estimated_Q11[1:n_sim]
)

# Tracer les probabilités de transition pour comparaison
# ggplot(comparison_Q11, aes(x = Time)) +
#   geom_line(aes(y = Simulated, color = "Simulated")) +
#   geom_line(aes(y = Estimated, color = "Estimated")) +
#   labs(title = "Comparison of Simulated and Estimated Q11(t)",
#        y = "Q11(t)",
#        color = "Legend") +
#   theme_minimal()
#   ggplot(comparison_Q22, aes(x = Time)) +
#   geom_line(aes(y = Simulated, color = "Simulated")) +
#   geom_line(aes(y = Estimated, color = "Estimated")) +
#   labs(title = "Comparison of Simulated and Estimated Q22(t)",
#        y = "Q22(t)",
#        color = "Legend") +
#   theme_minimal()
# Tracer les probabilités de transition pour comparaison
p_Q11 <- ggplot(comparison_Q11, aes(x = Time)) +
  geom_line(aes(y = Simulated, color = "Simulated")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "Comparison of Simulated and Estimated Q11(t)",
       y = "Q11(t)",
       color = "Legend") +
  theme_minimal()
ggsave("Q11_comparison.pdf", plot = p_Q11, width = 7, height = 5)

p_Q22 <- ggplot(comparison_Q22, aes(x = Time)) +
  geom_line(aes(y = Simulated, color = "Simulated")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "Comparison of Simulated and Estimated Q22(t)",
       y = "Q22(t)",
       color = "Legend") +
  theme_minimal()
ggsave("Q22_comparison.pdf", plot = p_Q22, width = 7, height = 5)



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
# ggplot(comparison_poisson[1:100,], aes(x = Time)) +
#   geom_line(aes(y = Simulated, color = "Simulated"), size = 1) +
#   geom_line(aes(y = Fitted, color = "Fitted"), size = 1, linetype = "dashed") +
#   labs(title = "Comparison of Simulated and Fitted Observations (Poisson)",
#        subtitle = "First 1000 Time Points",
#        y = "Count",
#        color = "Legend") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5),
#         legend.position = "bottom") +
#   scale_color_manual(values = c("Simulated" = "blue", "Fitted" = "red"))

# # Graphique pour les observations normales
# ggplot(comparison_normal[1:1000,], aes(x = Time)) +
#   geom_line(aes(y = Simulated, color = "Simulated"), size = 1) +
#   geom_line(aes(y = Fitted, color = "Fitted"), size = 1, linetype = "dashed") +
#   labs(title = "Comparison of Simulated and Fitted Observations (Normal)",
#        subtitle = "First 1000 Time Points",
#        y = "Value",
#        color = "Legend") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5),
#         legend.position = "bottom") +
#   scale_color_manual(values = c("Simulated" = "blue", "Fitted" = "red")
# Graphique pour les observations Poisson
p_poisson <- ggplot(comparison_poisson[1:100,], aes(x = Time)) +
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
ggsave("Poisson_comparison.pdf", plot = p_poisson, width = 7, height = 5)

# Graphique pour les observations normales
p_normal <- ggplot(comparison_normal[1:1000,], aes(x = Time)) +
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
  scale_color_manual(values = c("Simulated" = "blue", "Fitted" = "red"))
ggsave("Normal_comparison.pdf", plot = p_normal, width = 7, height = 5)

p3_bis <- ggplot(donnees, aes(x = date, y = temp, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Température extreme par état (données originales)", 
           x = "Time", y = "Température (°C)", color = "État") +
      theme_minimal()