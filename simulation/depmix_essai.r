# Charger la bibliothèque nécessaire
library(depmixS4)

# Définir les paramètres du modèle
set.seed(123)
n <- 365 * 3  # Nombre de jours simulés (3 ans de données)
states <- 2    # Nombre d'états
degree_obs_pol <- 1  # Degré des covariables trigonométriques
degree_trans_pol <- 1  # Degré des covariables trigonométriques
period <- 365  # Période saisonnière
time = 1:n

# Fonction pour calculer la matrice de transition Q avec dépendance saisonnière
compute_transition_matrix <- function(t, period, beta_matrix, K, degree) {
  P <- matrix(0, nrow = K, ncol = K)
  
  for (i in 1:K) {
    for (j in 2:K) {  # Éviter la diagonale principale
      P[i, j] <- beta_matrix[[i, j - 1]][1]  # Constante
      for (d in 1:degree) {  # Ajouter les termes trigonométriques
        P[i, j] <- P[i, j] + 
          beta_matrix[[i, j - 1]][2 * d] * cos(2 * pi * d * t / period) + 
          beta_matrix[[i, j - 1]][2 * d + 1] * sin(2 * pi * d * t / period)
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
simulate_states <- function(n, period, beta_matrix, K, degree) {
  states <- numeric(n)
  states[1] <- sample(1:n, 1)  # État initial
  
  for (t in 2:n) {
    Q_t <- compute_transition_matrix(t, period, beta_matrix, K, degree)
    states[t] <- sample(1:K, size = 1, prob = Q_t[states[t - 1], ])
  }
  return(states)
}

# Paramètres pour la transition saisonnière
beta <- list(
  c(1, 0.6, -0.4, 4, -4),     # Transitions de l'état 1 vers 2 d'autres
  c(0.4, 0.3, -0.3, 3, -3),   # Transitions de l'état 1 vers 3 d'autres
  c(-0.1, 0.1, 0.8, 5, -5),   # Transitions de l'état 2 vers 2 d'autres
  c(0.2, 0.3, -0.2, 2, -2),   # Transitions de l'état 2  vers 3 d'autres
  c(-0.5, 0.3, 0.9, 6, -6),  # Transitions de l'état 3 vers 2 d'autres,   
  c(0.4, -0.5, -0.1, 7, -7)  # Transitions de l'état 3 vers 3 d'autres
)

# Créer une matrice de coefficients pour les transitions
beta_matrix <- matrix(beta, nrow = states, ncol = states - 1, byrow = TRUE)

# Simuler les états
simulated_states <- simulate_states(n, period, beta_matrix, states, degree_trans_pol)

# Résumé des résultats
table(simulated_states)


# Simuler les observations avec \( m_k(t) \)
simulate_observations <- function(states, n, period, mu, delta, sd, degree) {
  observations <- numeric(n)
  for (t in 1:n) {
    k <- states[t]
    mean_t <- mu[k]
    for (d in 1:degree) {
      mean_t <- mean_t + delta[k, 2*d-1] * cos(2 * pi * d * t / period) + delta[k, 2*d] * sin(2 * pi * d * t / period)
    }
    observations[t] <- rnorm(1, mean = mean_t, sd = sd[k])
    
  }
  return(observations)
}

# Generalized function to compute trigonometric covariates
generate_trig_covariates <- function(time, period, degree) {
  trig_covs <- data.frame(time = time)
  for (d in 1:degree) {
    trig_covs[[paste0("sin_", d)]] <- sin(2 * pi * d * time / period)
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * time / period)
  }
  return(trig_covs)
}

# Paramètres pour les émissions
mu <- c(5, 10, -20)  # Moyennes des états
delta <- matrix(c(5, 4, 3, 2, -3, -2, -1, 0.1, 3, 8, 6, 2), nrow = states, byrow = TRUE)  # Coefficients pour cos et sin
sd <- c(1, 0.5, 5)  # Écarts-types

# Ajuster un modèle SHMM avec DepMixS4
trig_covs <- generate_trig_covariates(time, period, degree_trans_pol)


# Simulated response data
observations <- simulate_observations(simulated_states, n, period, mu, delta, sd, degree_obs_pol)
data <- cbind(
  obs = observations,
  trig_covs
)

# Define transition formula dynamically
transition_formula <- as.formula(
  paste("~", paste(names(trig_covs)[-1], collapse = " + "))
 )
#transition_formulas <- lapply(1:states, function(i) {
#  as.formula(paste("~", paste(names(trig_covs)[-1], collapse = " + ")))
# })
# Define obs formula dynamically
obs_formula <- as.formula(
  paste("obs ~", paste(names(trig_covs)[-1], collapse = " + "))
)


mod <- depmix(
  response = obs_formula,
  data = data,
  nstates = 2,
  family = gaussian(),
  transition = transition_formula,
  nstart = 10000)

fitted_model <- fit(mod, verbose = TRUE)

# Résultats
summary(fitted_model, which='transition')
summary(fitted_model, which='response')

simulate_observations_from_fitted <- function(fitted_model, states, n, period, degree) {
  observations <- numeric(n)
  for (t in 1:n) {
    k <- states[t]
    response_model <- fitted_model@response[[k]][[1]]@parameters
    response_params <- response_model$coefficients
    response_params_sd <- response_model$sd
    mean_t <- response_params[1]  # Intercept
    for (d in 1:degree) {
      mean_t <- mean_t + 
        response_params[2 * d] * cos(2 * pi * d * t / period) + 
        response_params[2 * d + 1] * sin(2 * pi * d * t / period)
    }
    observations[t] <- rnorm(1, mean = mean_t, sd = response_params_sd)
  }
  return(observations)
}

# Fonction pour simuler les états à partir des transitions estimées
simulate_states_from_fitted <- function(fitted_model, n, period, degree) {
  K <- fitted_model@nstates
  states <- numeric(n)
  states[1] <- 1  # Initial state
  
  # Generate trigonometric covariates
  trig_covs <- generate_trig_covariates(1:n, period, degree)
  
  for (t in 2:n) {
    # Extract transition matrix parameters
    transition_params <- summary(fitted_model, which = 'transition')
    
    # Compute transition probabilities
    P <- matrix(0, nrow = K, ncol = K)
    
    for (i in 1:K) {
      for (j in 1:K) {
        if (j == 1) {
          # Base probability (intercept)
          P[i, j] <- transition_params$coef[i, 1]
        } else {
          # Add trigonometric covariates
          for (d in 1:degree) {
            cos_coef <- transition_params$coef[i, 2*d]
            sin_coef <- transition_params$coef[i, 2*d+1]
            P[i, j] <- P[i, j] + 
              cos_coef * cos(2 * pi * d * t / period) + 
              sin_coef * sin(2 * pi * d * t / period)
          }
        }
      }
    }
    
    # Convert to transition probabilities using softmax
    Q <- exp(P)
    Q <- Q / rowSums(Q)
    
    # Sample next state
    states[t] <- sample(1:K, 1, prob = Q[states[t-1], ])
  }
  
  return(states)
}

# Simulate states using fitted model parameters
fitted_states <- simulate_states_from_fitted(
  fitted_model, 
  n, 
  period, 
  degree_trans_pol
)

# Simuler les états à partir des transitions estimées
fitted_states <- simulate_states_from_fitted(fitted_model, n, period, degree_trans_pol)

# Simuler les observations à partir des paramètres estimés
fitted_observations <- simulate_observations_from_fitted(fitted_model, simulated_states, n, period, degree_obs_pol)
# Comparer les observations générées par simulation et celles générées par les paramètres estimés
comparison <- data.frame(
  Time = 1:n,
  Simulated = observations,
  Fitted = fitted_observations
)

# Afficher les premières lignes de la comparaison
head(comparison)

# Tracer les observations pour comparaison
library(ggplot2)
ggplot(comparison, aes(x = Time)) +
  geom_line(aes(y = Simulated, color = "Simulated")) +
  geom_line(aes(y = Fitted, color = "Fitted")) +
  labs(title = "Comparison of Simulated and Fitted Observations",
       y = "Observations",
       color = "Legend") +
  theme_minimal()