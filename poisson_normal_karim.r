# Libraries
library(depmixS4)
library(tidyverse)
library(ggplot2)
library(readxl)
library(gridExtra)
library(stats)
# Définition de la nouvelle classe S4 pour customPoisson
setClass("customPoisson", 
         contains = "response",
         slots = c(
           period = "numeric"  # On ajoute un slot pour la période
         ))

# Définir une méthode générique pour customPoisson
setGeneric("customPoisson", 
           function(formula = response ~ trend + cos_1 + sin_1, 
                    data = NULL, 
                    pstart = NULL, 
                    fixed = NULL,
                    period = 52,
                    exposure_var = "exposure", 
                    ...) 
             standardGeneric("customPoisson"))

# Définir la méthode qui crée la classe de réponse
setMethod("customPoisson",
          signature(formula = "ANY"),
          function(formula = response ~ trend + cos_1 + sin_1, 
                   data = NULL, 
                   pstart = NULL,  
                   fixed = NULL,
                   period = 52,
                   exposure_var = "exposure", 
                   ...) {
            
            # Extraire les données du modèle
            if (!is.null(data)) {
              mf <- model.frame(formula, data = data)
              y <- as.matrix(model.response(mf))
              x <- model.matrix(formula, data = data)
              if (exposure_var %in% names(data)) {
                exposure <- data[[exposure_var]]
              } else {
                exposure <- rep(1, length(y))
                warning("Variable d'exposition non trouvée, utilisation de 1")
              }
            } else {
              y <- NULL
              x <- NULL
              exposure <- NULL
            }
            
            # Paramètres
            parameters <- list()
            npar <- 4  # lambda, alpha, gamma_cos, gamma_sin
            
            # Gestion des paramètres initiaux
            if (is.null(fixed)) fixed <- as.logical(rep(0, npar))
            if (!is.null(pstart)) {
              if (length(pstart) != npar) stop("length of 'pstart' must be ", npar)
              parameters$lambda <- pstart[1]      # Intercepte
              parameters$alpha <- pstart[2]       # Coefficient de tendance
              parameters$gamma_cos <- pstart[3]   # Coefficient cosinus
              parameters$gamma_sin <- pstart[4]   # Coefficient sinus
            } else {
              # Valeurs par défaut
              parameters$lambda <- -6
              parameters$alpha <- 0.01
              parameters$gamma_cos <- 0.2
              parameters$gamma_sin <- 0.1
            }
            
            # Créer l'objet customPoisson
            mod <- new("customPoisson",
                       parameters = parameters,
                       fixed = fixed,
                       x = x,
                       y = y,
                       npar = npar,
                       period = period)
            
            # Ajouter l'exposition comme attribut
            attr(mod, "exposure") <- exposure
            
            return(mod)
          })

# Méthode pour afficher l'objet
setMethod("show", "customPoisson",
          function(object) {
            cat("Model of type customPoisson (Poisson with seasonality and trend)\n")
            cat("Parameters:\n")
            cat("lambda (intercept): ", object@parameters$lambda, "\n")
            cat("alpha (trend): ", object@parameters$alpha, "\n")
            cat("gamma_cos (cosine): ", object@parameters$gamma_cos, "\n")
            cat("gamma_sin (sine): ", object@parameters$gamma_sin, "\n")
            cat("Period: ", object@period, "\n")
          })

# Méthode pour la densité
setMethod("dens", "customPoisson",
          function(object, log = FALSE) {
            # Extraction des paramètres
            lambda <- object@parameters$lambda
            alpha <- object@parameters$alpha
            gamma_cos <- object@parameters$gamma_cos
            gamma_sin <- object@parameters$gamma_sin
            
            # Extraction des covariables de la matrice de design
            if (!is.null(object@x)) {
              # On suppose que x contient l'intercepte et les covariables
              intercept <- object@x[, 1]
              
              # Vérifier si les colonnes existent avant de les extraire
              if (ncol(object@x) >= 2) trend <- object@x[, 2] else trend <- 0
              if (ncol(object@x) >= 3) cos_term <- object@x[, 3] else cos_term <- 0
              if (ncol(object@x) >= 4) sin_term <- object@x[, 4] else sin_term <- 0
              
              # Calcul de log(mu)
              log_mu <- lambda * intercept + 
                        alpha * trend + 
                        gamma_cos * cos_term + 
                        gamma_sin * sin_term
              
              mu <- exp(log_mu)
              
              # Appliquer l'exposition
              exposure <- attr(object, "exposure")
              if (!is.null(exposure)) {
                mu <- mu * exposure
              }
              
              # Calculer la log-vraisemblance Poisson
              dens_val <- dpois(object@y, lambda = mu, log = log)
              return(dens_val)
            } else {
              warning("Matrice de design non disponible")
              return(NA)
            }
          })

# Méthode pour obtenir les paramètres
setMethod("getpars", "customPoisson",
          function(object, which = "pars", ...) {
            switch(which,
                   "pars" = {
                     parameters <- numeric()
                     parameters <- unlist(object@parameters)
                     pars <- parameters
                   },
                   "fixed" = {
                     pars <- object@fixed
                   }
            )
            return(pars)
          })

# Méthode pour définir les paramètres
setMethod("setpars", "customPoisson",
          function(object, values, which = "pars", ...) {
            npar <- npar(object)
            if (length(values) != npar) stop("length of 'values' must be ", npar)
            
            # Déterminer si on définit les paramètres ou les contraintes fixes
            nms <- names(object@parameters)
            switch(which,
                   "pars" = {
                     object@parameters$lambda <- values[1]
                     object@parameters$alpha <- values[2]
                     object@parameters$gamma_cos <- values[3]
                     object@parameters$gamma_sin <- values[4]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters) <- nms
            return(object)
          })

# Méthode pour ajuster le modèle
setMethod("fit", "customPoisson",
          function(object, w) {
            if (missing(w)) w <- NULL
            
            # Extraction des données
            y <- object@y
            x <- object@x
            exposure <- attr(object, "exposure")
            
            if (is.null(exposure)) {
              exposure <- rep(1, length(y))
            }
            
            # Créer un dataframe pour l'ajustement
            fit_data <- data.frame(
              y = y,
              offset_val = log(exposure)
            )
            
            # Ajouter les covariables
            if (ncol(x) >= 2) fit_data$trend <- x[, 2]
            if (ncol(x) >= 3) fit_data$cos_term <- x[, 3]
            if (ncol(x) >= 4) fit_data$sin_term <- x[, 4]
            
            # Formule pour glm
            formula_str <- "y ~ "
            if (ncol(x) >= 2) formula_str <- paste0(formula_str, "trend + ")
            if (ncol(x) >= 3) formula_str <- paste0(formula_str, "cos_term + ")
            if (ncol(x) >= 4) formula_str <- paste0(formula_str, "sin_term + ")
            formula_str <- paste0(formula_str, "offset(offset_val)")
            
            # Ajuster le modèle GLM
            fit <- glm(as.formula(formula_str), 
                       data = fit_data, 
                       family = poisson(link = "log"),
                       weights = w)
            
            # Extraire les coefficients
            coefs <- coef(fit)
            
            # Mettre à jour les paramètres de l'objet
            pars <- numeric(4)
            pars[1] <- coefs[1]  # lambda (intercept)
            if (ncol(x) >= 2 && "trend" %in% names(coefs)) pars[2] <- coefs["trend"] else pars[2] <- 0
            if (ncol(x) >= 3 && "cos_term" %in% names(coefs)) pars[3] <- coefs["cos_term"] else pars[3] <- 0
            if (ncol(x) >= 4 && "sin_term" %in% names(coefs)) pars[4] <- coefs["sin_term"] else pars[4] <- 0
            
            # Mettre à jour l'objet
            object <- setpars(object, pars)
            return(object)
          })

# Méthode pour prédire
setMethod("predict", "customPoisson",
          function(object) {
            # Extraction des paramètres
            lambda <- object@parameters$lambda
            alpha <- object@parameters$alpha
            gamma_cos <- object@parameters$gamma_cos
            gamma_sin <- object@parameters$gamma_sin
            
            # Extraction des covariables
            if (!is.null(object@x)) {
              intercept <- object@x[, 1]
              
              # Vérifier si les colonnes existent avant de les extraire
              if (ncol(object@x) >= 2) trend <- object@x[, 2] else trend <- 0
              if (ncol(object@x) >= 3) cos_term <- object@x[, 3] else cos_term <- 0
              if (ncol(object@x) >= 4) sin_term <- object@x[, 4] else sin_term <- 0
              
              # Calcul de log(mu)
              log_mu <- lambda * intercept + 
                        alpha * trend + 
                        gamma_cos * cos_term + 
                        gamma_sin * sin_term
              
              mu <- exp(log_mu)
              
              # Appliquer l'exposition
              exposure <- attr(object, "exposure")
              if (!is.null(exposure)) {
                mu <- mu * exposure
              }
              
              return(mu)
            } else {
              warning("Matrice de design non disponible")
              return(NA)
            }
          })

setMethod("npar", "formula", 
          function(object) {
            # For formulas, return the number of parameters based on terms
            # This is just an example - adjust as needed
            rhs <- terms(object)
            return(length(attr(rhs, "term.labels")) + 1) # +1 for intercept
          })

# Méthode pour générer des valeurs aléatoires
setMethod("simulate", "customPoisson",
          function(object, nsim = 1, seed = NULL, times = 1) {
            if (!is.null(seed)) set.seed(seed)
            
            # Obtenir les valeurs prédites (mu)
            mu <- predict(object)
            
            # Générer des nombres aléatoires selon une distribution Poisson
            if (length(mu) > 0) {
              sim_values <- matrix(rpois(n = length(mu) * nsim, lambda = mu), 
                                   nrow = length(mu), 
                                   ncol = nsim)
              return(sim_values)
            } else {
              warning("Prédictions non disponibles")
              return(NULL)
            }
          })


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

# Exemple d'utilisation pour construire un modèle HMM avec customPoisson
# Exemple d'utilisation pour construire un modèle HMM avec customPoisson
build_hmm_with_custom_poisson <- function(data, K = 2, period = 52) {
  # Liste pour stocker les composantes de réponse pour chaque état
  resp_list <- list()

  simple_poisson <- glm(obs_poisson ~ trend + sin_1 + cos_1 + offset(log_exposure), 
                        family = poisson, data = data)
  simple_normal <- lm(obs_normal ~ trend + sin_1 + cos_1, data = data)
  
  # Extract coefficients as starting points
  poisson_coefs <- coef(simple_poisson)
  normal_coefs <- coef(simple_normal)
  normal_sd <- sd(residuals(simple_normal))
  
  # Modèle pour la composante Poisson (mortalité)
  mod_mort <- list()
  for (k in 1:K) {
    pstart_poisson <- c(
      lambda = poisson_coefs[1] * (1 + 0.1 * (k-1)),
      alpha = poisson_coefs[2] * (1 + 0.1 * (k-1)),
      gamma_cos = poisson_coefs[3] * (1 + 0.1 * (k-1)),
      gamma_sin = poisson_coefs[4] * (1 + 0.1 * (k-1))
    )
    
    mod_mort[[k]] <- customPoisson(
      formula = obs_poisson ~ trend + cos_1 + sin_1 ,
      data = data,
      pstart = pstart_poisson,
      period = period,
      exposure_var = "log_exposure"
    )
  }
  
  # Modèle pour la composante température (approche alternative)
  mod_temp <- list()
  for (k in 1:K) {
    pstart_temp <- c(
      intercept = normal_coefs[1] * (1 + 0.1 * (k-1)),
      trend = normal_coefs[2] * (1 + 0.1 * (k-1)),
      sin = normal_coefs[3] * (1 + 0.1 * (k-1)),
      cos = normal_coefs[4] * (1 + 0.1 * (k-1)),
      sd = normal_sd * (1 + 0.1 * (k-1))
    )
    mod_temp[[k]] <- GLMresponse(
      obs_normal ~ trend + sin_1 + cos_1,
      data = data,
      family = gaussian(),
      pstart = pstart
    )
  }
  
  # Construire la liste complète des réponses
  for (k in 1:K) {
    state_resp <- list()
    state_resp[[1]] <- mod_mort[[k]]
    state_resp[[2]] <- mod_temp[[k]]
    resp_list[[k]] <- state_resp
  }
  
  # Paramètres de transition initiaux
  trstart <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow = K)
  
  # Probabilités initiales
  instart <- rep(1/K, K)
  
  # Formule de transition
  tr_formula <- as.formula(paste("~", paste("cos_1", "sin_1", sep = " + ")))
  
  # Créer le modèle HMM
  mod <- makeDepmix(
    response = resp_list,          # Liste de listes de réponses par état
    transition = list(tr_formula), # Formule de transition
    prior = ~1,                    # Probabilité initiale
    homogeneous = FALSE,           # Transitions non-homogènes
    data = data,
    nstates = K,
    trstart = trstart,
    instart = instart
  )
  
  return(mod)
}

custom_mod <- build_hmm_with_custom_poisson(df, K = states)
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
