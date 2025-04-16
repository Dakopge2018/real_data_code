# Paramètres initiaux
set.seed(123)
period <- 365
T_total <- period * 5
K <- 2  # Nombre d'états
E <- 100  # Exposition constante

# Vrais paramètres
lambda_k <- c(-1, 0.5)  # Coefficients de base par état
alpha_k1 <- c(0.1, 0.2)  # Coefficients de tendance par état

# Matrice de transition (Q) et probabilités initiales (ξ)
# Ces valeurs seront fixées pour la simulation
Q_true <- matrix(c(0.95, 0.05, 0.1, 0.9), nrow = K, byrow = TRUE)
xi_true <- c(0.6, 0.4)

# Variables temporelles
time <- 1:T_total
a_t <- time %/% period

# Fonction pour générer une séquence d'états selon un processus de Markov
generate_states <- function(T_total, Q, xi) {
  states <- numeric(T_total)
  states[1] <- sample(1:nrow(Q), 1, prob = xi)
  
  for (t in 2:T_total) {
    states[t] <- sample(1:nrow(Q), 1, prob = Q[states[t-1], ])
  }
  
  return(states)
}

# Génération de la séquence d'états
states <- generate_states(T_total, Q_true, xi_true)

# Fonction pour générer les observations selon la séquence d'états
generate_observations <- function(states, lambda_k, alpha_k1, time, period, E, a_t) {
  T_total <- length(states)
  y <- numeric(T_total)
  
  # Générer des coefficients gamma pour chaque état
  gamma_k1 <- runif(K, -1, 1)
  gamma_k2 <- runif(K, -1, 1)
  
  # Stocker les paramètres gamma pour référence
  gammas <- list(gamma_k1 = gamma_k1, gamma_k2 = gamma_k2)
  
  for (t in 1:T_total) {
    k <- states[t]
    log_mu_k_t <- lambda_k[k] + 
                  gamma_k1[k] * cos(2 * pi * time[t] / period) + 
                  gamma_k2[k] * sin(2 * pi * time[t] / period)
    mu_k_t <- exp(log_mu_k_t)
    T_k_t <- alpha_k1[k] * a_t[t] * time[t]
    
    # Générer l'observation
    y[t] <- rpois(1, E * mu_k_t) + T_k_t
  }
  
  return(list(y = y, gammas = gammas))
}

# Génération des observations
observations <- generate_observations(states, lambda_k, alpha_k1, time, period, E, a_t)
y <- observations$y
gammas <- observations$gammas

# Créer le dataframe pour l'analyse
data <- data.frame(
  t = time,
  E = E,
  a_t = a_t,
  y = y,
  true_state = states,
  cos_t = cos(2 * pi * time / period),
  sin_t = sin(2 * pi * time / period)
)

library(dplyr)

# Fonction pour calculer la log-vraisemblance négative
neg_log_likelihood <- function(params, data, K, fixed_Q = NULL, fixed_xi = NULL) {
  # Extraction des paramètres
  lambda <- params[1:K]
  gamma1 <- params[(K+1):(2*K)]
  gamma2 <- params[(2*K+1):(3*K)]
  alpha <- params[(3*K+1):(4*K)]
  
  # Utiliser les matrices de transition et probabilités initiales fournies ou estimées
  if (is.null(fixed_Q)) {
    # Extraire les paramètres de la matrice de transition (log-odds)
    Q_params <- params[(4*K+1):(4*K+K*(K-1))]
    
    # Construire la matrice de transition Q
    Q <- matrix(0, nrow = K, ncol = K)
    idx <- 1
    for (i in 1:K) {
      for (j in 1:K) {
        if (i != j) {
          Q[i, j] <- Q_params[idx]
          idx <- idx + 1
        }
      }
    }
    
    # Appliquer la transformation softmax pour chaque ligne
    for (i in 1:K) {
      Q_exp <- exp(c(0, Q[i, -i]))  # 0 pour la diagonale (référence)
      Q_sum <- sum(Q_exp)
      Q[i, ] <- Q_exp / Q_sum
    }
  } else {
    Q <- fixed_Q
  }
  
  if (is.null(fixed_xi)) {
    # Extraire les probabilités initiales (log-odds)
    xi_params <- params[(4*K+K*(K-1)+1):(4*K+K*(K-1)+K-1)]
    
    # Calculer les probabilités initiales
    xi_exp <- exp(c(0, xi_params))  # 0 pour le premier état (référence)
    xi <- xi_exp / sum(xi_exp)
  } else {
    xi <- fixed_xi
  }
  
  # Calcul des log-probabilités d'émission pour chaque observation et chaque état
  log_emission <- matrix(0, nrow = nrow(data), ncol = K)
  
  for (k in 1:K) {
    log_mu <- lambda[k] + gamma1[k] * data$cos_t + gamma2[k] * data$sin_t
    mu <- exp(log_mu)
    trend <- alpha[k] * data$a_t * data$t
    
    # Valeur prédite totale pour chaque observation
    predicted <- data$E * mu + trend
    
    # Log-vraisemblance approximative pour chaque observation
    residuals <- data$y - predicted
    log_emission[, k] <- -residuals^2 / (2 * predicted) - 0.5 * log(2 * pi * predicted)
  }
  
  # Implémentation de l'algorithme Forward (en log-espace)
  T_total <- nrow(data)
  alpha_forward <- matrix(0, nrow = T_total, ncol = K)
  
  # Initialisation
  for (i in 1:K) {
    alpha_forward[1, i] <- log(xi[i]) + log_emission[1, i]
  }
  
  # Récursion
  for (t in 2:T_total) {
    for (j in 1:K) {
      max_val <- -Inf
      for (i in 1:K) {
        val <- alpha_forward[t-1, i] + log(Q[i, j])
        if (val > max_val) max_val <- val
      }
      
      sum_exp <- 0
      for (i in 1:K) {
        sum_exp <- sum_exp + exp(alpha_forward[t-1, i] + log(Q[i, j]) - max_val)
      }
      
      alpha_forward[t, j] <- max_val + log(sum_exp) + log_emission[t, j]
    }
  }
  
  # Log-vraisemblance du modèle
  log_likelihood <- max(alpha_forward[T_total, ])
  max_idx <- which.max(alpha_forward[T_total, ])
  
  for (i in 1:K) {
    if (i != max_idx) {
      log_likelihood <- log_likelihood + log(1 + exp(alpha_forward[T_total, i] - alpha_forward[T_total, max_idx]))
    }
  }
  
  return(-log_likelihood)  # Négatif car on minimise
}

# Algorithme de Viterbi
viterbi_algorithm <- function(params, data, K, fixed_Q = NULL, fixed_xi = NULL) {
  # Extraction des paramètres
  lambda <- params[1:K]
  gamma1 <- params[(K+1):(2*K)]
  gamma2 <- params[(2*K+1):(3*K)]
  alpha <- params[(3*K+1):(4*K)]
  
  # Utiliser les matrices de transition et probabilités initiales fournies ou estimées
  if (is.null(fixed_Q)) {
    # Extraire les paramètres de la matrice de transition (log-odds)
    Q_params <- params[(4*K+1):(4*K+K*(K-1))]
    
    # Construire la matrice de transition Q
    Q <- matrix(0, nrow = K, ncol = K)
    idx <- 1
    for (i in 1:K) {
      for (j in 1:K) {
        if (i != j) {
          Q[i, j] <- Q_params[idx]
          idx <- idx + 1
        }
      }
    }
    
    # Appliquer la transformation softmax pour chaque ligne
    for (i in 1:K) {
      Q_exp <- exp(c(0, Q[i, -i]))  # 0 pour la diagonale (référence)
      Q_sum <- sum(Q_exp)
      Q[i, ] <- Q_exp / Q_sum
    }
  } else {
    Q <- fixed_Q
  }
  
  if (is.null(fixed_xi)) {
    # Extraire les probabilités initiales (log-odds)
    xi_params <- params[(4*K+K*(K-1)+1):(4*K+K*(K-1)+K-1)]
    
    # Calculer les probabilités initiales
    xi_exp <- exp(c(0, xi_params))  # 0 pour le premier état (référence)
    xi <- xi_exp / sum(xi_exp)
  } else {
    xi <- fixed_xi
  }
  
  # Calcul des log-probabilités d'émission pour chaque observation et chaque état
  log_emission <- matrix(0, nrow = nrow(data), ncol = K)
  
  # Calculer les valeurs prédites pour chaque état
  predicted_values <- array(0, dim = c(nrow(data), K, 2))  # 2 composantes: Poisson et Tendance
  
  for (k in 1:K) {
    log_mu <- lambda[k] + gamma1[k] * data$cos_t + gamma2[k] * data$sin_t
    mu <- exp(log_mu)
    trend <- alpha[k] * data$a_t * data$t
    
    # Valeur prédite totale pour chaque observation
    predicted <- data$E * mu + trend
    
    # Stocker les composantes prédites
    predicted_values[, k, 1] <- data$E * mu  # Composante Poisson
    predicted_values[, k, 2] <- trend        # Composante Tendance
       # Log-vraisemblance approximative pour chaque observation
    r0..0000.
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    5--esiduals <- data$y - predicted
    log_emission[, k] <- -residuals^2 / (2 * predicted) - 0.5 * log(2 * pi * predicted)
  }
  
  # Début de l'algorithme de Viterbi
  T_total <- nrow(data)
  delta <- matrix(0, nrow = T_total, ncol = K)
  psi <- matrix(0, nrow = T_total, ncol = K)
  0
  # Initialisation (t=1)
  for (i in 1:K) {
    delta[1, i] <- log(xi[i]) + log_emission[1, i]
    psi[1, i] <- 0  # Pas de prédécesseur à t=1
  }
  
  # Récursion (t=2,...,T)
  for (t in 2:T_total) {
    for (j in 1:K) {
      max_val <- -Inf
      max_idx <- 0
      
      for (i in 1:K) {
        val <- delta[t-1, i] + log(Q[i, j])
        if (val > max_val) {
          max_val <- val
          max_idx <- i
        }
      }
      
      delta[t, j] <- max_val + log_emission[t, j]
      psi[t, j] <- max_idx
    }
  }
  
  # Termination
  best_path_prob <- max(delta[T_total, ])
  best_path_last_state <- which.max(delta[T_total, ])
  
  # Path backtracking
  best_path <- numeric(T_total)
  best_path[T_total] <- best_path_last_state
  
  for (t in (T_total-1):1) {
    best_path[t] <- psi[t+1, best_path[t+1]]
  }
  
  # Calculer les prédictions pour le meilleur chemin
  predictions <- data.frame(
    t = data$t,
    y_actual = data$y,
    y_pred = numeric(T_total),
    poisson_comp = numeric(T_total),
    trend_comp = numeric(T_total),
    estimated_state = best_path,
    true_state = data$true_state
  )
  
  for (t in 1:T_total) {
    k <- best_path[t]
    predictions$y_pred[t] <- predicted_values[t, k, 1] + predicted_values[t, k, 2]
    predictions$poisson_comp[t] <- predicted_values[t, k, 1]
    predictions$trend_comp[t] <- predicted_values[t, k, 2]
  }
  
  return(list(
    predictions = predictions,
    best_path = best_path,
    best_path_prob = best_path_prob,
    Q = Q,
    xi = xi
  ))
}

# Pour la démonstration, nous allons fixer Q et xi aux vraies valeurs
# Dans une application réelle, on pourrait les estimer également
fixed_Q <- Q_true
fixed_xi <- xi_true

# Valeurs initiales des paramètres (seulement pour lambda, gamma1, gamma2, alpha)
initial_params <- c(
  lambda = c(-0.5, 0.5),    # K valeurs
  gamma1 = c(0, 0),         # K valeurs
  gamma2 = c(0, 0),         # K valeurs
  alpha = c(0.05, 0.15)     # K valeurs
)

# Définition des contraintes: A %*% theta + B >= 0
# Contraintes:
# 1. alpha_1 > 0
# 2. alpha_2 > 0
# 3. lambda_2 > lambda_1

A <- matrix(0, nrow = 3, ncol = length(initial_params))
A[1, 3*K+1] <- 1      # alpha_1 > 0
A[2, 3*K+2] <- 1      # alpha_2 > 0
A[3, 2] <- 1          # lambda_2 > 0
A[3, 1] <- -1         # lambda_1 < 0

B <- c(0, 0, 0)

# Vérifier que le point initial satisfait les contraintes
if (any(A %*% initial_params + B < 0)) {
  stop("Le point initial ne satisfait pas les contraintes!")
}

# Optimisation avec contraintes
opt_result <- constrOptim(
  theta = initial_params,
  f = function(params) neg_log_likelihood(params, data, K, fixed_Q, fixed_xi),
  grad = NULL,  # Pas de gradient analytique
  ui = A,
  ci = B,
  control = list(maxit = 10000)
)

# Affichage des résultats
estimated_params <- opt_result$par
param_names <- c(paste0("lambda_", 1:K), 
                 paste0("gamma1_", 1:K), 
                 paste0("gamma2_", 1:K), 
                 paste0("alpha_", 1:K))
names(estimated_params) <- param_names

# Obtenir les vrais paramètres
true_params <- c(
  lambda_k,
  gammas$gamma_k1,
  gammas$gamma_k2,
  alpha_k1
)
names(true_params) <- param_names

# Affichage des résultats
results_df <- data.frame(
  Parameter = names(estimated_params),
  Estimated = estimated_params,
  True = true_params,
  Difference = estimated_params - true_params,
  Percent_Error = abs((estimated_params - true_params) / true_params * 100)
)

print(results_df)
cat("\nConvergence:", ifelse(opt_result$convergence == 0, "Successful", "Failed"), "\n")
cat("Log-likelihood:", -opt_result$value, "\n")
cat("\nContraintes respectées:", all(A %*% estimated_params + B >= 0), "\n")

# Exécuter l'algorithme de Viterbi avec les paramètres estimés
viterbi_results <- viterbi_algorithm(estimated_params, data, K, fixed_Q, fixed_xi)

# Évaluer la précision de la classification des états
classification_accuracy <- mean(viterbi_results$best_path == data$true_state)
cat("\nPrécision de la classification des états:", classification_accuracy, "\n")

# Calculer le R² pour les prédictions
predictions <- viterbi_results$predictions
r2 <- 1 - sum((predictions$y_actual - predictions$y_pred)^2) / 
            sum((predictions$y_actual - mean(predictions$y_actual))^2)
cat("R² global:", r2, "\n")

# Visualisation des résultats
library(ggplot2)

# Graphique des observations vs prédictions
p1 <- ggplot(predictions, aes(x = t)) +
  geom_line(aes(y = y_actual), alpha = 0.4, color = "gray") +
  geom_line(aes(y = y_pred, color = factor(estimated_state)), size = 1) +
  labs(title = "Observations vs Prédictions",
       subtitle = paste("R² =", round(r2, 4), "| Précision états =", round(classification_accuracy, 4)),
       x = "Temps (jours)", y = "Valeur",
       color = "État estimé") +
  theme_minimal()

# Comparaison des états estimés vs réels
p2 <- ggplot(predictions, aes(x = t)) +
  geom_line(aes(y = true_state, color = "Réel"), size = 1) +
  geom_line(aes(y = estimated_state, color = "Estimé"), size = 1, linetype = "dashed") +
  scale_y_continuous(breaks = 1:K) +
  labs(title = "États cachés: réels vs estimés",
       x = "Temps (jours)", y = "État",
       color = "Type") +
  theme_minimal()

# Décomposition des composantes pour les états estimés
p3 <- ggplot(predictions, aes(x = t, color = factor(estimated_state))) +
  geom_line(aes(y = poisson_comp), size = 1) +
  labs(title = "Composante Poisson estimée",
       x = "Temps (jours)", y = "Valeur",
       color = "État estimé") +
  theme_minimal()

p4 <- ggplot(predictions, aes(x = t, color = factor(estimated_state))) +
  geom_line(aes(y = trend_comp), size = 1) +
  labs(title = "Composante Tendance estimée",
       x = "Temps (jours)", y = "Valeur",
       color = "État estimé") +
  theme_minimal()

print(p1)
print(p2)
print(p3)
print(p4)