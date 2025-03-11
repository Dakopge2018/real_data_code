# Charger la bibliothèque nécessaire
library(depmixS4)
library(ggplot2)

# Définir les paramètres du modèle
set.seed(123)
n <- 365 * 5  # Nombre de jours simulés (3 ans de données)
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

# Paramètres pour la transition saisonnière
beta <- list(
  c(1, 0.7, 0.5),     # Transitions de l'état 1 vers 2 d'autres
  c(-1, -0.6, 0.7)   # Transitions de l'état 2 vers 2 d'autres
)

# Créer une matrice de coefficients pour les transitions
beta_matrix <- matrix(beta, nrow = states, ncol = states - 1, byrow = TRUE)

# Simuler les états
simulated_states <- f_simulate_states(n, period, beta_matrix, states, degree_trans_pol)

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
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * time / period)}
  return(trig_covs)
}

# Paramètres pour les émissions
mu <- c(-1, 2)  # Moyennes des états
delta <- matrix(c(2.5, 4 , -1.5, 3.5), nrow = states, byrow = TRUE)  # Coefficients pour cos et sin
sd <- c(1, 0.25)  # Écarts-types

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
  transition = transition_formula)

set.seed(1)
fitted_model  <- multistart(mod,
  nstart = 50,  # 10 initialisations différentes
  initIters = 30,  # 10 itérations EM pour chaque initialisation
  emcontrol = em.control(
    maxit = 50000,  # Max 500 itérations EM
    tol = 1e-08,  # Tolérance pour convergence
    crit = "relative",  # Critère de convergence
    random.start = TRUE,  # Randomisation des initialisations
    classification = "Hard"  # Classification 
    ))


# Résultats
summary(fitted_model, which='transition')
summary(fitted_model, which='response')

simulate_observations_from_fitted <- function(fitted_model, states, n, period, degree) {
  observations <- numeric(n)
  for (t in 1:n) {
    k <- states[t]
    #if (k == 1) {
     # k <- 2
    #}
    #else if (k == 2) {
    #  k <- 1
   # }
    response_model <- fitted_model@response[[k]][[1]]@parameters
    response_params <- response_model$coefficients
    response_params_sd <- response_model$sd
    mean_t <- response_params[1]  # Intercept
    for (d in 1:degree) {
      mean_t <- mean_t + 
        response_params[2 * d+1] * cos(2 * pi * d * t / period) + 
        response_params[2 * d] * sin(2 * pi * d * t / period)
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

# Simuler les états à partir des transitions estimées
# fitted_states <- simulate_states_from_fitted(fitted_model, n, period, degree_trans_pol)
fitted_states = posterior(fitted_model)$state
# Simuler les observations à partir des paramètres estimés
fitted_observations <- simulate_observations_from_fitted(fitted_model, fitted_states, n, period, degree_obs_pol)
# Comparer les observations générées par simulation et celles générées par les paramètres estimés
comparison <- data.frame(
  Time = 1:n,
  Simulated = observations,
  Fitted = fitted_observations
)

# Afficher les premières lignes de la comparaison
head(comparison)

# Tracer les observations pour comparaison
ggplot(comparison[1:100,], aes(x = Time)) +
  geom_line(aes(y = Simulated, color = "Simulated")) +
  geom_line(aes(y = Fitted, color = "Fitted")) +
  labs(title = "Comparison of Simulated and Fitted Observations",
       y = "Observations",
       color = "Legend") +
  theme_minimal()

# Calculer des statistiques de comparaison
mse <- mean((comparison$Simulated - comparison$Fitted)^2)
mae <- mean(abs(comparison$Simulated - comparison$Fitted))
cor_coef <- cor(comparison$Simulated, comparison$Fitted)

# Afficher les métriques de comparaison
cat("Mean Squared Error:", mse, "\n")
cat("Mean Absolute Error:", mae, "\n")
cat("Correlation Coefficient:", cor_coef, "\n")

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

# Appliquer l'algorithme Forward-Backward
#fb_results <- forwardbackward(fitted_model)

# Récupérer la matrice des probabilités de transition à chaque instant t
#xi <- fb_results$xi  # Dimensions: (T x nstates x nstates)
#estimated_Q12 = xi[,1,2]
#estimated_Q21 = xi[,2,1]
#estimated_Q11 = xi[,1,1]
#estimated_Q22 = xi[,2,2]

# Extraire les probabilités de transition estimées
estimated_beta <- list(
    list(fitted_model@transition[[1]]@parameters$coefficients[,2]),
    list(fitted_model@transition[[2]]@parameters$coefficients[,2])
    )

est_beta_matrix <- matrix(estimated_beta, nrow = states, ncol = states - 1, byrow = TRUE)
# Comparer les probabilités de transition simulées et estimées
estimated_Q11 <- sapply(1:n, function(t) compute_simulated_transition_probs(t, period, est_beta_matrix, states, degree_trans_pol)[1, 1])
estimated_Q22 <- sapply(1:n, function(t) compute_simulated_transition_probs(t, period, est_beta_matrix, states, degree_trans_pol)[2, 2])
estimated_Q12 <- sapply(1:n, function(t) compute_simulated_transition_probs(t, period, est_beta_matrix, states, degree_trans_pol)[1, 2])
estimated_Q21 <- sapply(1:n, function(t) compute_simulated_transition_probs(t, period, est_beta_matrix, states, degree_trans_pol)[2, 1])
comparison_Q11 <- data.frame(
  Time = 1:700,
  Simulated = simulated_Q11[1:700],
  Estimated = estimated_Q11[1:700]
)

comparison_Q22 <- data.frame(
  Time = 1:700,
  Simulated = simulated_Q22[1:700],
  Estimated = estimated_Q22[1:700]
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

state1 <- rep(1, n)
state2 <- rep(2, n)

observations_state1 <- simulate_observations(state1, n, period, mu, delta, sd, degree_obs_pol)
observations_state2 <- simulate_observations(state2, n, period, mu, delta, sd, degree_obs_pol)

# Combine the data into a single data frame
realizations <- data.frame(
  Time = rep(1:n, 2),
  Observations = c(observations_state1, observations_state2),
  State = factor(rep(c("State 1", "State 2"), each = n))
)

# Plot the data using ggplot2
ggplot(realizations, aes(x = Time, y = Observations, color = State)) +
  geom_line() +
  labs(title = "Realization of Y1 to Y1000 for State 1 and State 2",
       x = "Time",
       y = "Observations",
       color = "State") +
  theme_minimal()