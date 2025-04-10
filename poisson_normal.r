# Libraries
library(depmixS4)
library(tidyverse)
library(glm)
library(ggplot2)
library(readxl)
library(gridExtra)
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
data <-na.omit(data)
data <- data[is.finite(rowSums(data[c('year','No_year', 'Temperature', 'Death_counts', 'Weekly_exposure')])), ]

# Vérifier s'il y a des NaN dans la colonne Death_counts
if (any(is.infinite(log(data$Weekly_exposure)))) {
  cat("Il y a des valeurs NaN dans la colonne Death_counts.\n")
} else {
  cat("Il n'y a pas de valeurs NaN dans la colonne Death_counts.\n")
}



generate_trig_covariates <- function(week_no, year_no, period, degree) {
  trig_covs <- data.frame(time = week_no, trend = year_no)
  for (d in 1:degree) {
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * week_no / period)
    trig_covs[[paste0("sin_", d)]] <- sin(2 * pi * d * week_no / period)}
  return(trig_covs)
}
trig_covs <- generate_trig_covariates(data$No_week, data$No_year,period, degree_trans_pol)

df = data.frame(
    obs_poisson = data$Death_counts,
    obs_normal = data$Temperature,
    rate = data$Death_counts / data$Weekly_exposure,
    log_exposure = log(data$Weekly_exposure),
    exposure = data$Weekly_exposure,
    trend = data$No_year,
    trig_covs)

# Calculer les percentiles 5% et 95% pour la colonne obs_poisson
quantiles_obs_poisson <- quantile(df$obs_poisson, probs = c(0.01, 0.99))

# Filtrer les données pour retirer les valeurs en dehors des percentiles 5% et 95% pour obs_poisson
df <- subset(df, obs_poisson >= quantiles_obs_poisson[1] & obs_poisson <= quantiles_obs_poisson[2])

# Supprimer les lignes avec des valeurs manquantes ou infinies
df <- na.omit(df)
df <- df[complete.cases(df), ]
df <- df[is.finite(rowSums(df)), ]

# visualisation des données
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
  labs(title="Exposure", x="Temps (semaine)", y="Exposure") +
  theme_minimal()


p4 <- ggplot(df, aes(x=obs_normal, y=obs_poisson)) +
  geom_point() +
  geom_smooth(method="loess") +
  labs(title="Relation entre température et décès", 
       x="Température (°C)", y="Nombre de décès") +
  theme_minimal()

grid.arrange(p1, p2, p3,p4, ncol=1)

# Transition formula
transition_formula <- as.formula(paste("~", paste(names(trig_covs)[c(-1,-2)], collapse = " + ")))

# Model
mod <- depmix(obs_poisson ~ offset(exposure) + sin_1 + cos_1 , 
  data = df,
  maxiter = 3,
  nstates = states,
  family = poisson(link = "log")
)
fitted_mod = fit(mod)
set.seed(1)

# fitted_model  <- multistart(mod,
#   nstarts = 20,  # 10 initialisations
#   initIters = 10, # 10 itérations EM pour chaque initialisation
#   emcontrol = em.control(
#     maxit = 500,  # Max 500 itérations EM
#     tol = 1e-08,  # Tolérance pour convergence
#     crit = "relative",  # Critère de convergence
#     random.start = TRUE,  # Randomisation des initialisations
#     classification = "hard"  # Classification 
#     ))
summary(fitted_mod, which = 'transition')
summary(fitted_mod, which = 'response')

# Récupération des états les plus probables
etats_predits <- posterior(fitted_mod, type = 'viterbi')
df$etat <- etats_predits$state

# Visualisation des données avec états prédits
p4 <- ggplot(df, aes(x=time, y=factor(etat), color=factor(etat))) +
  geom_point(size=3) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="États cachés prédits", x="Temps (semiane)", y="État", color="État") +
  theme_minimal()

p5 <- ggplot(df, aes(x=time, y=obs_poisson, color=factor(etat))) +
  geom_line() +
  geom_point(size=2) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="Décès par état", x="Temps (semaine)", y="Nombre de décès", color="État") +
  theme_minimal()

p6 <- ggplot(df, aes(x=time, y=obs_normal, color=factor(etat))) +
  geom_line() +
  geom_point(size=2) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="Température par état", x="Temps (semaine)", y="Température (°C)", color="État") +
  theme_minimal()

grid.arrange(p4, p5, p6, ncol=1)

model_poisson <- glm((obs_poisson ~offset(log_exposure) +trend + sin_1 + cos_1 ),
  data = df,
  family = poisson(link = "log"),
  nstart = 10,
  control = list(maxit = 1000))


summary(model_poisson)