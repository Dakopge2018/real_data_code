# Libraries
library(depmixS4)
library(ggplot2)

# Model parameters
set.seed(123)
n <- 365 * 10
states <- 2
degree_obs_pol <- 1
degree_trans_pol <- 1
period <- 365
time = 1:n

# Fonction pour calculer la matrice de transition Q avec dépendance saisonnière
compute_transition_matrix <- function(t, period, beta_matrix, K, degree) {
  P <- matrix(0, nrow = K, ncol = K)
  
  for (i in 1:K) {
    for (j in 2:K) {  # Éviter la diagonale principale
      P[i, j] <- beta_matrix[[i, j - 1]][1]  # Constante
      for (d in 1:degree) {  # Ajouter les termes trigonométriques
        P[i, j] <- P[i, j] + 
          beta_matrix[[i, j - 1]][2 * d] * sin(2 * pi * d * t / period) + 
          beta_matrix[[i, j - 1]][2 * d + 1] * cos(2 * pi * d * t / period)
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

# Fonction pour simuler les états
f_simulate_states <- function(n, period, beta_matrix, K, degree) {
  states <- numeric(n)
  states[1] <- sample(1:K, 1)  # État initial
  
  for (t in 2:n) {
    Q_t <- compute_transition_matrix(t, period, beta_matrix, K, degree)
    states[t] <- sample(1:K, size = 1, prob = Q_t[states[t - 1], ])
  }
  return(states)
}
# Generalized function to compute trigonometric covariates
generate_trig_covariates <- function(time, period, degree) {
  trig_covs <- data.frame(time = time)
  for (d in 1:degree) {
    cos_term <- cos(2 * pi * d * time / period)
    sin_term <- sin(2 * pi * d * time / period)
    trig_covs[[paste0("cos_", d)]] <- cos_term
    trig_covs[[paste0("sin_", d)]] <- sin_term
    # trig_covs[[paste0("Exp(cos_", d, "+sin_", d, ")")]] <- exp(cos_term + sin_term)
  }
  return(trig_covs)
}

# Modified observation simulation for Poisson
simulate_observations <- function(states, n, period, mu, delta, degree) {
 observations <- numeric(n)
 for (t in 1:n) {
   k <- states[t]
   lambda_t <- exp(mu[k])  # Base rate
   for (d in 1:degree) {
     lambda_t <- lambda_t * exp(delta[k, 2*d-1] * sin(2 * pi * d * t / period) + 
                               delta[k, 2*d] * cos(2 * pi * d * t / period))
   }
   observations[t] <- rpois(1, lambda = lambda_t)
 }
 return(observations)
}

# Parameters
beta <- list(
 c(1, 0.7, 0.5),
 c(-1, -0.6, 0.7)
)
beta_matrix <- matrix(beta, nrow = states, ncol = states - 1, byrow = TRUE)

# Modified emission parameters for Poisson
mu <- c(2, 0.7)  # Log baseline rates
delta <- matrix(c(1, 0.7, -0.9, 0.1), nrow = states, byrow = TRUE)

# Simulate data
simulated_states <- f_simulate_states(n, period, beta_matrix, states, degree_trans_pol)
trig_covs <- generate_trig_covariates(time, period, degree_trans_pol)
observations <- simulate_observations(simulated_states, n, period, mu, delta, degree_obs_pol)

# Prepare data
data <- cbind(
 obs = observations,
 trig_covs
)

# Fit model with Poisson response
transition_formula <- as.formula(paste("~", paste(names(trig_covs)[-1], collapse = " + ")))
obs_formula <- as.formula(paste("obs ~", paste(names(trig_covs)[-1], collapse = " + ")))

mod <- depmix(
 response = obs_formula,
 data = data,
 nstates = 2,
 family = poisson(),  # Changed to Poisson
 transition = transition_formula)

set.seed(1)
fitted_model  <- multistart(mod,
  nstart = 1000,  # 10 initialisations différentes
  initIters = 20,  # 10 itérations EM pour chaque initialisation
  emcontrol = em.control(
    maxit = 500,  # Max 500 itérations EM
    tol = 1e-08,  # Tolérance pour convergence
    crit = "relative",  # Critère de convergence
    random.start = TRUE,  # Randomisation des initialisations
    classification = "Hard"  # Classification 
    ))

fitted_model <- fit(mod, verbose = TRUE)
summary(fitted_model, which='transition')
summary(fitted_model, which='response')


# Modified simulation from fitted model for Poisson
simulate_observations_from_fitted <- function(fitted_model, states, n, period, degree) {
 observations <- numeric(n)
 for (t in 1:n) {
   k <- ifelse(states[t] == 1, 2, 1)
   response_model <- fitted_model@response[[k]][[1]]@parameters
   response_params <- response_model$coefficients
   
   lambda_t <- exp(response_params[1])  # Intercept
   for (d in 1:degree) {
     lambda_t <- lambda_t * exp(response_params[2 * d] * sin(2 * pi * d * t / period) + 
                               response_params[2 * d + 1] * cos(2 * pi * d * t / period))
   }
   observations[t] <- rpois(1, lambda = lambda_t)
 }
 return(observations)
}

# Rest of comparison and plotting code remains the same
# fitted_states <- simulate_states_from_fitted(fitted_model, n, period, degree_trans_pol)
fitted_states = posterior(fitted_model)$state
fitted_observations <- simulate_observations_from_fitted(fitted_model, fitted_states, n, period, degree_obs_pol)

comparison <- data.frame(
 Time = 1:n,
 States = simulated_states,
 Simulated = observations,
 Fitted = fitted_observations
)

# Plotting and metrics
ggplot(comparison[1:100,], aes(x = Time)) +
 geom_line(aes(y = Simulated, color = "Simulated")) +
 geom_line(aes(y = Fitted, color = "Fitted")) +
 labs(title = "Comparison of Simulated and Fitted Observations (Poisson)",
      y = "Count",
      color = "Legend") +
 theme_minimal()

# Calculate metrics
mse <- mean((comparison$Simulated - comparison$Fitted)^2)
mae <- mean(abs(comparison$Simulated - comparison$Fitted))
cor_coef <- cor(comparison$Simulated, comparison$Fitted)
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
  Simulated = simulated_Q11,
  Estimated = estimated_Q11
)

comparison_Q22 <- data.frame(
  Time = 1:200,
  Simulated = simulated_Q22,
  Estimated = estimated_Q22
)

# Tracer les probabilités de transition pour comparaison
ggplot(comparison_Q11, aes(x = Time)) +
  geom_line(aes(y = Simulated, color = "Simulated")) +
  geom_line(aes(y = Estimated, color = "Estimated")) +
  labs(title = "Comparison of Simulated and Estimated Q11(t)",
       y = "Q11(t)",
       color = "Legend") +
  theme_minimal()
