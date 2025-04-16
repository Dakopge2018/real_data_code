# Libraries
library(depmixS4)
library(tidyverse)
library(ggplot2)
library(readxl)
library(gridExtra)
library(stats)

# Définition d'une réponse Poisson personnalisée pour la mortalité
# où log(mu) = lambda + alpha*trend + gamma_cos*cos + gamma_sin*sin
customPoisson <- function(formula = response ~ trend + sin_1 + cos_1, data = NULL, pstart = NULL, 
                         fixed = NULL, exposure_var = "exposure", ...) {
  # Extraction des variables
  if (is.null(data)) {
    response <- NULL
    trend <- NULL
    sin_1 <- NULL
    cos_1 <- NULL
    exposure <- NULL
  } else {
    response <- model.response(model.frame(formula, data = data))
    covs <- model.matrix(formula, data = data)
    exposure <- data[[exposure_var]]
  }
  
  # Paramètres initiaux
  if (is.null(pstart)) {
    pstart <- c(lambda = -6, alpha = 0.01, gamma_cos = 0.2, gamma_sin = 0.1)
  }
  
  # Paramètres fixes
  if (is.null(fixed)) {
    fixed <- c(FALSE, FALSE, FALSE, FALSE)
  }
  
  # Nombre de paramètres
  npar <- length(pstart)
  
  # Fonction de densité
  dens <- function(x, pars) {
    # Extraction des paramètres
    lambda <- pars[1]
    alpha <- pars[2]
    gamma_cos <- pars[3]
    gamma_sin <- pars[4]
    
    # Matrice de design pour les covariables
    covs <- cbind(1, model.matrix(formula, data = data)[,-1]) # Intercepter + covariables
    
    # Calcul de log(mu) selon le modèle spécifié
    log_mu <- covs %*% c(lambda, alpha, gamma_cos, gamma_sin)
    mu <- exp(log_mu)
    
    # Application du facteur d'exposition
    mu <- exposure * mu
    
    # Calcul de la log-vraisemblance
    dens <- dpois(x, lambda = mu, log = TRUE)
    
    return(dens)
  }
  
  # Fonction de prédiction
  pred <- function(pars) {
    # Extraction des paramètres
    lambda <- pars[1]
    alpha <- pars[2]
    gamma_cos <- pars[3]
    gamma_sin <- pars[4]
    
    # Matrice de design pour les covariables
    covs <- cbind(1, model.matrix(formula, data = data)[,-1]) # Intercepter + covariables
    
    # Calcul de log(mu) selon le modèle spécifié
    log_mu <- covs %*% c(lambda, alpha, gamma_cos, gamma_sin)
    mu <- exp(log_mu)
    
    # Application du facteur d'exposition
    return(exposure * mu)
  }
  
  # Fonction de génération aléatoire
  rval <- function(pars) {
    # Utiliser la fonction de prédiction pour obtenir mu
    mu <- pred(pars)
    
    # Générer des valeurs aléatoires selon une distribution Poisson
    return(rpois(length(mu), lambda = mu))
  }
  
  # Construction de l'objet de réponse
  res <- list(
    npar = npar,
    dens = dens,
    pred = pred,
    rval = rval,
    pars = pstart,
    fixed = fixed,
    formula = formula,
    name = "customPoisson"
  )
  
  class(res) <- "response"
  return(res)
}

# Model parameters
set.seed(123)
states <- 2
degree_obs_pol <- 1
degree_trans_pol <- 1
period <- 52

# Read the Excel file
# Replace 'path_to_your_file.xlsx' with the actual path to your Excel file
data <- read_excel("main_database.xlsx", sheet = "database")

# Check for and handle missing or infinite values
data <- na.omit(data)
data <- data[is.finite(rowSums(data[c('year','No_year', 'Temperature', 'Death_counts', 'Weekly_exposure')])), ]

# Vérifier s'il y a des NaN dans la colonne Death_counts
if (any(is.infinite(log(data$Weekly_exposure)))) {
  cat("Il y a des valeurs NaN dans la colonne Death_counts.\n")
} else {
  cat("Il n'y a pas de valeurs NaN dans la colonne Death_counts.\n")
}

# Générer les covariables trigonométriques
generate_trig_covariates <- function(week_no, year_no, period, degree) {
  trig_covs <- data.frame(time = week_no, trend = year_no)
  for (d in 1:degree) {
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * week_no / period)
    trig_covs[[paste0("sin_", d)]] <- sin(2 * pi * d * week_no / period)}
  return(trig_covs)
}
trig_covs <- generate_trig_covariates(data$No_week, data$No_year, period, degree_trans_pol)

# Préparation des données
df = data.frame(
    obs_poisson = data$Death_counts,
    obs_normal = data$Temperature,
    rate = data$Death_counts / data$Weekly_exposure,
    log_exposure = log(data$Weekly_exposure),
    exposure = data$Weekly_exposure,
    trend = data$No_year,
    trig_covs)

# Calculer les percentiles 1% et 99% pour la colonne obs_poisson
#quantiles_obs_poisson <- quantile(df$obs_poisson, probs = c(0.01, 0.99))

# Filtrer les données pour retirer les valeurs en dehors des percentiles 1% et 99% pour obs_poisson
#df <- subset(df, obs_poisson >= quantiles_obs_poisson[1] & obs_poisson <= quantiles_obs_poisson[2])

# Supprimer les lignes avec des valeurs manquantes ou infinies
df <- na.omit(df)
df <- df[complete.cases(df), ]
df <- df[is.finite(rowSums(df)), ]

# Visualisation des données
p1 <- ggplot(df, aes(x=time, y=obs_poisson)) +
  geom_line() +
  geom_point() +
  labs(title="Décès observés", x="Temps (semaine)", y="Nombre de décès") +
  theme_minimal()

p2 <- ggplot(df, aes(x=time, y=obs_normal)) +
  geom_line() +
  geom_point() +
  labs(title="Température observée", x="Temps (semaine)", y="Température (°C)") +
  theme_minimal()

p3 <- ggplot(df, aes(x=time, y=exposure)) +
  geom_line() +
  geom_point() +
  labs(title="Exposition", x="Temps (semaine)", y="Exposition") +
  theme_minimal()

p4 <- ggplot(df, aes(x=obs_normal, y=obs_poisson)) +
  geom_point() +
  geom_smooth(method="loess") +
  labs(title="Relation entre température et décès", 
       x="Température (°C)", y="Nombre de décès") +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, ncol=2)

# Construction du modèle HMM avec customPoisson
build_custom_hmm <- function(data, K = 2) {
  # Modèle pour la composante mortalité (customPoisson)
  mod_mortality <- list()
  for (k in 1:K) {
    # Paramètres initiaux pour chaque état
    pstart <- c(
      lambda = -5.5 - 0.5 * (k-1),  # Différent pour chaque état
      alpha = 0.01 + 0.01 * (k-1),
      gamma_cos = 0.2 + 0.1 * (k-1),
      gamma_sin = 0.1 + 0.05 * (k-1)
    )
    
    mod_mortality[[k]] <- customPoisson(
      formula = obs_poisson ~ trend + sin_1 + cos_1,
      data = data,
      pstart = pstart,
      exposure_var = "exposure"
    )
  }
  
  # Modèle pour la composante température (Normale)
  mod_temp <- list()
  for (k in 1:K) {
    pstart <- c(
      intercept = 15 - 5 * (k - 1),     # Intercept différent par état
      trend = 0.05 - 0.02 * (k - 1),    # Coeff du trend
      sin = 8 + 2 * (k - 1),            # Saisonnalité sin
      cos = 2 + (k - 1),                # Saisonnalité cos
      sd = 2 + (k*0.5)                            # Ecart-type
    )
    mod_temp[[k]] <- GLMresponse(
      obs_normal ~ trend + sin_1 + cos_1,
      data = data,
      family = gaussian(),
      pstart = pstart
      )
    }
  
  # Paramètres de transition initiaux
  trstart <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow = K)
  
  # Initialisation des probabilités des états initiaux
  instart <- rep(1/K, K)
  
  # Formule de transition
  tr_formula <- as.formula(paste("~", paste("sin_1", "cos_1", sep = " + ")))
  
  # Construction du modèle HMM avec makedepmix
  mod <- makeDepmix(
    response = list(
      list(mod_mortality[[1]], mod_mortality[[2]]),  # Réponses Poisson pour chaque état
      list(mod_temp[[1]], mod_temp[[2]])          # Réponses normales pour chaque état
    ),
    transition = list(tr_formula),         # Transition dépendant des covariables saisonnières
    prior = ~1,                            # Probabilité initiale constante
    homogeneous = FALSE,                   # Transitions non-homogènes
    data = data,
    nstates = K,
    trstart = trstart,
    instart = instart
  )
  
  return(mod)
}

# Créer et ajuster le modèle HMM personnalisé
custom_mod <- build_custom_hmm(df, K = states)
fit_start_time <- Sys.time()
custom_fit <- fit(custom_mod, verbose = TRUE, emc = em.control(maxit = 500, tol = 1e-08))
fit_end_time <- Sys.time()
cat("Temps d'ajustement du modèle:", difftime(fit_end_time, fit_start_time, units = "mins"), "minutes\n")

# Résumé du modèle
summary(custom_fit)
summary(custom_fit, which = 'transition')
summary(custom_fit, which = 'response')

# Récupération des états les plus probables
etats_predits <- posterior(custom_fit, type = 'viterbi')
df$etat <- etats_predits$state

# Visualisation des données avec états prédits
p5 <- ggplot(df, aes(x=time, y=factor(etat), color=factor(etat))) +
  geom_point(size=3) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="États cachés prédits", x="Temps (semaine)", y="État", color="État") +
  theme_minimal()

p6 <- ggplot(df, aes(x=time, y=obs_poisson, color=factor(etat))) +
  geom_line() +
  geom_point(size=2) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="Décès par état", x="Temps (semaine)", y="Nombre de décès", color="État") +
  theme_minimal()

p7 <- ggplot(df, aes(x=time, y=obs_normal, color=factor(etat))) +
  geom_line() +
  geom_point(size=2) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="Température par état", x="Temps (semaine)", y="Température (°C)", color="État") +
  theme_minimal()

p8 <- ggplot(df, aes(x=obs_normal, y=obs_poisson, color=factor(etat))) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="Relation température-mortalité par état", 
       x="Température (°C)", y="Nombre de décès", color="État") +
  theme_minimal()

grid.arrange(p5, p6, p7, p8, ncol=2)

# Extraction des paramètres estimés
params <- getpars(custom_fit)
print("Paramètres estimés du modèle:")
print(params)

# Fonction pour calculer et visualiser les prédictions du modèle
calculate_predictions <- function(model_fit, data) {
  # Extraire les états prédits
  states_pred <- posterior(model_fit, type = 'viterbi')$state
  
  # Extraire les paramètres pour chaque composante et état
  params <- getpars(model_fit)
  
  # Préparer le dataframe pour les prédictions
  pred_df <- data
  pred_df$state <- states_pred
  pred_df$mort_pred <- NA
  pred_df$temp_pred <- NA
  
  # Pour chaque observation, calculer les prédictions selon l'état
  for (i in 1:nrow(pred_df)) {
    state_i <- pred_df$state[i]
    
    # Indice pour les paramètres de mortalité selon l'état
    if (state_i == 1) {
      # Paramètres pour état 1 (à adapter selon la structure exacte)
      mort_params <- params[c(1:4)]
    } else {
      # Paramètres pour état 2 (à adapter selon la structure exacte)
      mort_params <- params[c(5:8)]
    }
    
    # Calculer log(mu) pour la mortalité
    covs_mort <- c(1, pred_df$trend[i], pred_df$sin_1[i], pred_df$cos_1[i])
    log_mu <- sum(covs_mort * mort_params)
    mu <- exp(log_mu) * pred_df$exposure[i]
    
    pred_df$mort_pred[i] <- mu
    
    # De même pour la température (à adapter)
    # ...
  }
  
  # Visualiser les prédictions vs observations
  p_mort <- ggplot(pred_df, aes(x=time)) +
    geom_line(aes(y=obs_poisson), color="black") +
    geom_line(aes(y=mort_pred), color="red") +
    labs(title="Mortalité observée vs prédite", 
         x="Temps", y="Décès", 
         color="") +
    theme_minimal()
  
  print(p_mort)
  
  return(pred_df)
}

# Essayer de calculer les prédictions (peut nécessiter des ajustements)
# pred_df <- calculate_predictions(custom_fit, df)

# Sauvegarde du modèle et des résultats
saveRDS(custom_fit, "hmm_customPoisson_model.rds")
saveRDS(df, "hmm_results_data.rds")