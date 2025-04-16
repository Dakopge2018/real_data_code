# Libraries
library(depmixS4)
library(tidyverse)
library(ggplot2)
library(readxl)
library(gridExtra)
library(stats)
library(LaplacesDemon)

set.seed(123)
states <- 2
degree_obs_pol <- 1
degree_trans_pol <- 1
period <- 52
# Lecture des données depuis un fichier Excel
data <- read_excel("main_database.xlsx", sheet = "database")

df = data.frame(
    obs_poisson = data$Death_counts,
  #obs_poisson = round(data$Death_counts/100,0),
    obs_normal = data$Temperature,
    rate = data$Death_counts / data$Weekly_exposure,
    log_exposure = log(data$Weekly_exposure),
    exposure = data$Weekly_exposure,
    trend = scale(data$No_year))
# Prétraitement des données
df <- df %>%
  mutate(
    time = data$No_week,
    cos_1 = cos(2 * pi * time / 52.18),
    sin_1 = sin(2 * pi * time / 52.18),
  ) %>%
  drop_na() %>%
  filter_all(all_vars(is.finite(.)))

# Vérification des colonnes requises
required_cols <- c("obs_poisson", "obs_normal", "trend", "cos_1", "sin_1", "exposure")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) stop("Colonnes manquantes : ", paste(missing_cols, collapse = ", "))


# build_hmm_with_exposure <- function(data, K = 2, period = 52) {
#   # Liste pour stocker les composantes de réponse pour chaque état
#   resp_list <- list()
  
#   # Modèle pour la composante Poisson (mortalité) avec l'exposition comme offset
#   for (k in 1:K) {
#     pstart_poisson <- c(
#       lambda = 0.9 - 0.5 * (k-1),
#       alpha = 0.01 + 0.01 * (k-1),
#       gamma_cos = -0.2 + 0.1 * (k-1),
#       gamma_sin = 0.1 + 0.05 * (k-1)
#     )

#     obs_formula_death <- as.formula(paste("obs_poisson ~ trend + sin_1 + cos_1"))
#     obs_formula_death_rate <- as.formula(paste("rate ~ trend + sin_1 + cos_1"))
#     # Ajout de l'exposition comme offset dans le modèle Poisson
#     mod_mort <- GLMresponse(
#       #obs_formula_death,  # Formule avec covariables
#       obs_formula_death_rate,
#       data = data,
#       family = gaussian(),       # Nombre d'états
#       #family = poisson(link = "log"),
#       #pstart = pstart_poisson,
#       #offset = log(data$exposure)  # L'exposition est incluse comme offset
#     )
    
#     resp_list[[k]] <- mod_mort
#   }
  
#   # Modèle pour la composante température (approche alternative)
#   for (k in 1:K) {
#     pstart_temp <- c(
#       intercept = 15 - 5 * (k - 1),     # Intercept différent par état
#       trend = 0.05 - 0.02 * (k - 1),    # Coeff du trend
#       sin = 8 + 2 * (k - 1),            # Saisonnalité sin
#       cos = 2 + (k - 1),                # Saisonnalité cos
#       sd = 2 + (k*0.5)                  # Ecart-type
#     )
#     obs_formula_temp <- as.formula(paste("obs_normal ~ trend + sin_1 + cos_1"))


#     mod_temp <- GLMresponse(
#       obs_formula_temp,
#       data = data,
#       family = gaussian(),
#       pstart = pstart_temp
#     )
    
#     resp_list[[k]] <- c(resp_list[[k]], mod_temp)
#   }

#   # Paramètres de transition initiaux
#   # trstart <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow = K)
  
#   # Probabilités initiales
#   # instart <- rep(1/K, K)
#   # Adjust transition probabilities
#   pstart1 <- c(logit(0.6), 0, 0, logit(0.4), 0, 0)
#   pstart2 <- c(logit(0.4), 0, 0, logit(0.6), 0, 0)

#   transition <- list()
#   transition[[1]] <- transInit(~cos_1+sin_1,nstates=2,data=data, pstart=pstart1)
#   transition[[2]] <- transInit(~cos_1+sin_1,nstates=2,data=data, pstart=pstart2)

#   instart=c(0.1,0.9)
#   inMod <- transInit(~ 1, ns = 2, ps = instart,data = data.frame(1))
  
#   # Formule de transition (simplifiée)
#   tr_formula <- ~ as.formula(cos_1 + sin_1) # Formule de transition avec covariables trigonométriques
  
#   # Créer le modèle HMM
#   mod <- makeDepmix(
#     response = resp_list,          # Liste de modèles pour chaque état
#     transition = transition, # Formule de transition
#     #prior = inMod,                    # Probabilité initiale
#     homogeneous = FALSE,           # Transitions non-homogènes
#     data = data,
#     nstates = K
#     #trstart = trstart,
#     #instart = instart
#   )
  
#   return(mod)
# }
#simple_poisson <- glm(obs_poisson ~ trend + sin_1 + cos_1 + offset(log_exposure), 
#                family = poisson(link ='log'), data = df)
                

build_hmm_with_exposure <- function(data, K = 2, period = 52) {
  # Liste pour stocker les composantes de réponse pour chaque état
  resp_list <- list()
  
  # Fit simple GLMs to get reasonable starting values
  simple_poisson <- glm(obs_poisson ~ trend + sin_1 + cos_1 + offset(log_exposure), family = poisson(link ='log'), data = data) # nolint
  simple_normal <- lm(obs_normal ~ trend + sin_1 + cos_1, data = data)
  
  # Extract coefficients as starting points
  poisson_coefs <- coef(simple_poisson)
  normal_coefs <- coef(simple_normal)
  normal_sd <- sd(residuals(simple_normal))
  
  # Modèle pour la composante Poisson (mortalité) avec l'exposition comme offset
  for (k in 1:K) {
    # Slightly different starting values for each state
    pstart_poisson <- c(
      lambda = poisson_coefs[1] * (1 + 0.1 * (k-1)),
      alpha = poisson_coefs[2] * (1 + 0.1 * (k-1)),
      gamma_cos = poisson_coefs[3] * (1 + 0.1 * (k-1)),
      gamma_sin = poisson_coefs[4] * (1 + 0.1 * (k-1))
    )

    obs_formula_death <- as.formula("obs_poisson ~ offset(log(exposure)) + trend + sin_1 + cos_1")
    # Ajout de l'exposition comme offset dans le modèle Poisson
    mod_mort <- GLMresponse(
      obs_formula_death,
      data = data,
      family = poisson(link = "log"),
      pstart = pstart_poisson
     # offset = data$log_exposure  # Use pre-computed log_exposure
    )
    
    resp_list[[k]] <- mod_mort
  }
  
  # Modèle pour la composante température (approche alternative)
  for (k in 1:K) {
    pstart_temp <- c(
      intercept = normal_coefs[1] * (1 + 0.1 * (k-1)),
      trend = normal_coefs[2] * (1 + 0.1 * (k-1)),
      sin = normal_coefs[3] * (1 + 0.1 * (k-1)),
      cos = normal_coefs[4] * (1 + 0.1 * (k-1)),
      sd = normal_sd * (1 + 0.1 * (k-1))
    )
    
    obs_formula_temp <- as.formula("obs_normal ~ trend + sin_1 + cos_1")

    mod_temp <- GLMresponse(
      obs_formula_temp,
      data = data,
      family = gaussian(),
      pstart = pstart_temp
    )
    
    resp_list[[k]] <- c(resp_list[[k]], mod_temp)
  }

  # Simpler transition model initialization
  transition <- list()
  for (k in 1:K) {
    # Less extreme starting values for transition probabilities
    trProbs <- rep(0.2, K)
    trProbs[k] <- 0.8  # Higher probability of staying in the same state
    
    pstart_trans <- c(qlogis(trProbs[-k]), rep(0, 5))
    
    transition[[k]] <- transInit(~cos_1+sin_1, nstates=K, data=data, pstart=pstart_trans)
  }

  # Simple initial state probabilities
  instart <- rep(1/K, K)
  inMod <- transInit(~ 1, ns = K, ps = instart, data = data.frame(1))
  
  # Créer le modèle HMM
  mod <- makeDepmix(
    response = resp_list,
    transition = transition,
    prior = inMod,
    homogeneous = FALSE,
    data = data,
    nstates = K
  )
  
  return(mod)
}
# First try a simpler model with fewer iterations
fit_start_time <- Sys.time()
custom_mod <- build_hmm_with_exposure(df, K = states)
tryCatch({
  custom_fit <- fit(custom_mod, verbose = TRUE, 
                   emc = em.control(maxit = 100, tol = 1e-06))
  
  # If successful, print summary
  summary(custom_fit)
}, error = function(e) {
  cat("Error in fitting:", e$message, "\n")
  # Suggest further simplification if needed
  # cat("Consider further simplifying the model or checking data\n")
})
fit_end_time <- Sys.time()
cat("Fitting time:", difftime(fit_end_time, fit_start_time, units = "mins"), "minutes\n")

# custom_mod <- build_hmm_with_exposure(df, K = states)
# fit_start_time <- Sys.time()
# custom_fit <- fit(custom_mod, verbose = TRUE, emc = em.control(maxit = 500, tol = 1e-08))
fit_end_time <- Sys.time()
cat("Temps d'ajustement du modèle:", difftime(fit_end_time, fit_start_time, units = "mins"), "minutes\n")

# Résumé du modèle
summary(custom_fit)
summary(custom_fit, which = 'transition')
summary(custom_fit, which = 'response')