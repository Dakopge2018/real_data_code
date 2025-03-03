# Libraries
library(depmixS4)
library(ggplot2)
library(readxl)

# Model parameters
set.seed(123)
states <- 2
degree_obs_pol <- 1
degree_trans_pol <- 1
period <- 52

# Read the Excel file
data <- read_excel("main_database.xlsx", sheet = "database")

# Check for and handle missing or infinite values
data <- na.omit(data)
data <- data[is.finite(rowSums(data[c('year','No_year', 'Temperature', 'Death_counts', 'Weekly_exposure')])), ]

# Fonction pour générer les covariables trigonométriques
generate_trig_covariates <- function(week_no, year_no, period, degree) {
  trig_covs <- data.frame(time = week_no, trend = year_no)
  for (d in 1:degree) {
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * week_no / period)
    trig_covs[[paste0("sin_", d)]] <- sin(2 * pi * d * week_no / period)
  }
  return(trig_covs)
}

# Génération des covariables trigonométriques
trig_covs <- generate_trig_covariates(data$No_week, data$No_year, period, degree_obs_pol)

# Préparation du dataframe pour le modèle
df <- data.frame(
    obs_poisson = data$Death_counts,
    obs_normal = data$Temperature,
    exposure = data$Weekly_exposure,
    trend = data$No_year,
    trig_covs)

# Vérification et nettoyage des données
df <- na.omit(df)
df <- df[is.finite(rowSums(df)), ]

# Vérification des valeurs extrêmes
summary(df)

# Formule pour la transition entre états
# Correction: utiliser les noms appropriés des colonnes trigonométriques
transition_formula <- as.formula(paste("~", paste(c("sin_1", "cos_1"), collapse = " + ")))

# S'assurer que la formule est bien construite
print(transition_formula)

# Définition du modèle
mod <- depmix(
  list(
    obs_poisson ~ trend + sin_1 + cos_1 + offset(log(exposure)), 
    obs_normal ~ sin_1 + cos_1 + trend
  ),
  data = df,
  nstates = states,
  family = list(poisson(), gaussian()),
  transition = transition_formula
)

# Au lieu d'utiliser makeDepmixConstraints, voici comment définir des contraintes:
# Par exemple, pour fixer certains paramètres de transition
# 1. Obtenir les paramètres libres actuels
freepars <- getpars(mod)
# 2. Créer un vecteur d'égal longueur avec TRUE/FALSE pour indiquer quels paramètres sont libres
# Par exemple, pour fixer les paramètres de transition du premier état:
fixed <- rep(TRUE, length(freepars))
# Les indices des paramètres dépendent de votre modèle spécifique
# Nous devons donc les identifier
# Exemple: supposons que les paramètres de transition commencent à l'indice 20
transition_start_idx <- 20  # Ceci est hypothétique! Ajuster selon votre modèle
fixed[transition_start_idx:(transition_start_idx+5)] <- FALSE  # Exemple

# Définir des contraintes directement lors de l'ajustement
try({
  # Première tentative avec optimisation standard
  fitted_mod <- fit(mod, verbose = TRUE, 
                   emcontrol = em.control(maxit = 500, tol = 1e-8))
})

# Si la première tentative échoue, essayer avec une meilleure initialisation
if(!exists("fitted_mod") || fitted_mod@message != "convergence") {
  message("Premier ajustement échoué, tentative avec initialisation aléatoire...")
  # Réinitialiser des valeurs aléatoires pour les paramètres initiaux
  set.seed(456)  # Différente graine
  # Méthode alternative: utiliser donlp2 qui peut être plus robuste
  fitted_mod <- fit(mod, verbose = TRUE, 
                   emcontrol = em.control(maxit = 1000, tol = 1e-10),
                   method = "donlp2")
}

# Vérifier la convergence
if(fitted_mod@message != "convergence") {
  warning("Le modèle n'a pas convergé correctement. Essayez différentes initialisations ou paramètres.")
}

# Extraire et afficher les paramètres
parameters <- getpars(fitted_mod)
print(parameters)

# Extraire les états cachés estimés
states_df <- data.frame(
  Time = 1:length(fitted_mod@posterior$state),
  State = fitted_mod@posterior$state,
  Probability = apply(fitted_mod@posterior$S, 1, max)
)

# Visualisation des états cachés
ggplot(states_df, aes(x = Time, y = State, color = factor(State))) + 
  geom_point() + 
  theme_minimal() + 
  ggtitle("Estimated Hidden States Over Time") +
  scale_color_discrete(name = "État")

# Évaluation du modèle
print(paste("BIC:", BIC(fitted_mod)))
print(paste("AIC:", AIC(fitted_mod)))

# Sauvegarde du modèle
saveRDS(fitted_mod, "fitted_hmm_model.rds")