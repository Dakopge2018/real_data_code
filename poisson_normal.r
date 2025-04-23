# Bibliothèques nécessaires
library(hmmTMB)
library(tidyverse)
library(ggplot2)
library(readxl)
library(gridExtra)
library(stats)
pal <- hmmTMB:::hmmTMB_cols

set.seed(123)
n_states <- 2
period <- 52.18
d <- 2
# Lecture des données depuis un fichier Excel
df <- read_excel("main_database.xlsx", sheet = "database")

# Préparation des données
data <- data.frame(
  death_counts = df$Death_counts,
  temperature = df$Temperature,
  exposure = df$Weekly_exposure,
  trend = df$No_year,
  t = df$No_week
) %>%
  mutate(
    cos_1 = cos(2 * pi * t / 52.18),
    sin_1 = sin(2 * pi * t / 52.18),
    log_exposure = log(exposure)
  ) %>%
  drop_na() %>%
  filter_all(all_vars(is.finite(.)))

# Vérifier la structure des données
str(data)
# Créer les composantes saisonnières
  # Nombre de termes dans le polynôme trigonométrique
season <- matrix(0, nrow = nrow(data), ncol = 2 * d)
for (l in 1:d) {
  season[, 2 * l - 1] <- cos(2 * pi * l * data$t / period)
  season[, 2 * l] <- sin(2 * pi * l * data$t / period)
}

data <- cbind(data, season)
colnames(data)[(ncol(data) - 2 * d + 1):ncol(data)] <- paste0("season_", 1:(2 * d))
# Processus d'état caché
hid <- MarkovChain$new(data = df, n_states = 2)

# Ajuster des modèles simples pour obtenir des estimations initiales
simple_poisson <- glm(death_counts ~ trend + season_1 + season_2 + season_3 + season_4 + offset(log(exposure)), 
                      family = poisson(link = 'log'), data = data)
simple_normal <- lm(temperature ~ trend + season_1 + season_2 + season_3 + season_4, data = data)

# Extraire les coefficients des modèles simples
coef_poisson <- coef(simple_poisson)
coef_normal <- coef(simple_normal)

# Extraire les valeurs ajustées pour créer des estimations différenciées par état
fitted_mortality <- fitted(simple_poisson)
fitted_temperature <- fitted(simple_normal)

# Calculer des moyennes séparées pour les estimations initiales des deux états
mortality_cutoff <- median(fitted_mortality)
mortality_low <- mean(fitted_mortality[fitted_mortality <= mortality_cutoff])
mortality_high <- mean(fitted_mortality[fitted_mortality > mortality_cutoff])

temp_cutoff <- median(fitted_temperature)
temp_low <- mean(fitted_temperature[fitted_temperature <= temp_cutoff])
temp_high <- mean(fitted_temperature[fitted_temperature > temp_cutoff])

temp_sd_low <- sd(data$temperature[fitted_temperature <= temp_cutoff])
temp_sd_high <- sd(data$temperature[fitted_temperature > temp_cutoff])

# Utiliser ces estimations dans le modèle HMM
obs <- Observation$new(
  data = data,
  dists = list(death_counts = "pois", temperature = "norm"),
  par = list(
    death_counts = list(
      rate = c(mortality_low, mortality_high),
      beta = list(
        "(Intercept)" = c(coef_poisson["(Intercept)"], coef_poisson["(Intercept)"]),
        "trend" = c(coef_poisson["trend"], coef_poisson["trend"]),
        #"trend" = c(coef_poisson["trend"], coef_poisson["trend"])
        "season_1" = c(coef_poisson["season_1"], coef_poisson["season_1"]),
        "season_2" = c(coef_poisson["season_2"], coef_poisson["season_2"]),
        "season_3" = c(coef_poisson["season_3"], coef_poisson["season_3"]),
        "season_4" = c(coef_poisson["season_4"], coef_poisson["season_4"])
      )
    ),
    temperature = list(
      mean = c(temp_low, temp_high),
      sd = c(temp_sd_low, temp_sd_high),
      beta = list(
        "(Intercept)" = c(coef_normal["(Intercept)"], coef_normal["(Intercept)"]),
        "trend" = c(coef_normal["trend"], coef_normal["trend"]),
        "season_1" = c(coef_normal["season_1"], coef_normal["season_1"]),
        "season_2" = c(coef_normal["season_2"], coef_normal["season_2"]),
        "season_3" = c(coef_normal["season_3"], coef_normal["season_3"]),
        "season_4" = c(coef_normal["season_4"], coef_normal["season_4"])
      )
    )
  ),
  formulas = list(
    death_counts = list(rate = ~ trend + season_1 + season_2 + season_3 + season_4 + offset(log(exposure))),
    #death_counts = list(mean = ~ trend + offset(log(exposure))),
    temperature = list(mean = ~ trend + season_1 + season_2 + season_3 + season_4)
  ),
  n_states = n_states
)
# Combiner les modèles
hmm <- HMM$new(obs = obs, hid = hid)

# Ajuster le modèle avec plus d'options de convergence
hmm$fit(silent = TRUE, control = list(maxit = 1000))

# Résumé du modèle
print(hmm)
summary(hmm)

# Tracer les probabilités d'état
hmm$plot_ts()

# Visualiser les distributions d'observation
hmm$plot_dist("death_counts")
hmm$plot_dist("temperature")

# Résidus pseudo pour les deux variables
pr <- hmm$pseudores()
par(mfrow = c(2, 2))
qqnorm(pr$mortality, main = "QQ-plot - Mortalité")
abline(0, 1)
acf(pr$mortality, main = "ACF - Mortalité")
qqnorm(pr$temperature, main = "QQ-plot - Température")
abline(0, 1)
acf(pr$temperature, main = "ACF - Température")