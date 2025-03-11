# Chargement des packages nécessaires
library(depmixS4)
library(ggplot2)
library(gridExtra)

#------------------------------------------------
# PARTIE 1: PRÉPARATION ET EXPLORATION DES DONNÉES
#------------------------------------------------

# Supposons que nous avons des données réelles de décès et température mensuelles
# Chargement des données (à remplacer par vos propres données)
# données_exemple <- read.csv("chemin/vers/vos/donnees.csv")

# Création d'un exemple de données (à remplacer par vos données réelles)
set.seed(123)
n <- 36  # 3 ans de données mensuelles
temps <- 1:n
annees <- rep(2021:2023, each=12)
mois <- rep(1:12, 3)

# Création de covariables saisonnières
cos_t <- cos(2*pi*temps/12)
sin_t <- sin(2*pi*temps/12)

# Simulation de données pour l'exemple (à remplacer par vos données réelles)
# Température avec cycle saisonnier
temp <- 15 + 10*cos_t + rnorm(n, 0, 2)
# Décès avec effet saisonnier et dépendance à la température
deces <- rpois(n, lambda = exp(3 + 0.5*cos_t - 0.02*temp))

# Création du dataframe
donnees_reelles <- data.frame(
  temps = temps,
  annee = annees,
  mois = mois,
  cos_t = cos_t,
  sin_t = sin_t,
  deces = deces,
  temperature = temp
)

# Exploration des données
head(donnees_reelles)
summary(donnees_reelles)

# Visualisation des données réelles
p1 <- ggplot(donnees_reelles, aes(x=temps, y=deces)) +
  geom_line() +
  geom_point() +
  labs(title="Décès observés", x="Temps (mois)", y="Nombre de décès") +
  theme_minimal()

p2 <- ggplot(donnees_reelles, aes(x=temps, y=temperature)) +
  geom_line() +
  geom_point() +
  labs(title="Température observée", x="Temps (mois)", y="Température (°C)") +
  theme_minimal()

p3 <- ggplot(donnees_reelles, aes(x=temperature, y=deces)) +
  geom_point() +
  geom_smooth(method="loess") +
  labs(title="Relation entre température et décès", 
       x="Température (°C)", y="Nombre de décès") +
  theme_minimal()

grid.arrange(p1, p2, p3, ncol=1)

#------------------------------------------------
# PARTIE 2: DÉFINITION ET AJUSTEMENT DU MODÈLE HMM
#------------------------------------------------

# Nombre d'états à définir (ici nous supposons 2 états: été/hiver)
nb_etats <- 2

# Définition du modèle HMM
# - Variable 1: décès (distribution de Poisson)
# - Variable 2: température (distribution Gaussienne)
# - Les transitions entre états dépendent des covariables saisonnières

modele_hmm <- depmix(list(deces ~ 1, temperature ~ 1), 
                    data=donnees_reelles,
                    nstates=nb_etats,
                    family=list(poisson(), gaussian()),
                    transition=~cos_t+sin_t)

# Visualisation de la structure du modèle avant ajustement
summary(modele_hmm)

# Ajustement du modèle (estimation des paramètres)
# Note: l'option 'emcontrol=em.control(maxit=500)' augmente le nombre d'itérations
# pour assurer la convergence
cat("Ajustement du modèle en cours, veuillez patienter...\n")
modele_ajuste <- fit(modele_hmm, emcontrol=em.control(maxit=500, tol=1e-8))

# Vérification de la convergence
cat("Convergence du modèle:", ifelse(modele_ajuste@message=="Log likelihood converged to within tol. (relative change)", 
                          "Réussie", paste("Échouée:", modele_ajuste@message)), "\n")

# Affichage des résultats de l'estimation
summary(modele_ajuste)

# Récupération des états les plus probables
etats_predits <- posterior(modele_ajuste)
donnees_reelles$etat <- etats_predits$state

# Visualisation des données avec états prédits
p4 <- ggplot(donnees_reelles, aes(x=temps, y=factor(etat), color=factor(etat))) +
  geom_point(size=3) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="États cachés prédits", x="Temps (mois)", y="État", color="État") +
  theme_minimal()

p5 <- ggplot(donnees_reelles, aes(x=temps, y=deces, color=factor(etat))) +
  geom_line() +
  geom_point(size=2) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="Décès par état", x="Temps (mois)", y="Nombre de décès", color="État") +
  theme_minimal()

p6 <- ggplot(donnees_reelles, aes(x=temps, y=temperature, color=factor(etat))) +
  geom_line() +
  geom_point(size=2) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="Température par état", x="Temps (mois)", y="Température (°C)", color="État") +
  theme_minimal()

grid.arrange(p4, p5, p6, ncol=1)

# Calcul des probabilités de transition conditionnelles
getpars(modele_ajuste)[grep("transition", names(getpars(modele_ajuste)))]

#------------------------------------------------
# PARTIE 3: SIMULATION À PARTIR DU MODÈLE AJUSTÉ
#------------------------------------------------

# Simulation de nouvelles données à partir du modèle ajusté
set.seed(456)
donnees_simulees <- simulate(modele_ajuste, nsim=1)

# Préparation pour visualisation
donnees_sim <- data.frame(
  temps = 1:n,
  mois = rep(1:12, 3),
  annee = rep(2024:2026, each=12),  # Simulation pour les années suivantes
  deces = donnees_simulees$deces,
  temperature = donnees_simulees$temperature,
  etat = donnees_simulees$state,
  cos_t = cos_t,
  sin_t = sin_t
)

# Visualisation des données simulées
p7 <- ggplot(donnees_sim, aes(x=temps, y=factor(etat), color=factor(etat))) +
  geom_point(size=3) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="États cachés simulés", x="Temps (mois)", y="État", color="État") +
  theme_minimal()

p8 <- ggplot(donnees_sim, aes(x=temps, y=deces, color=factor(etat))) +
  geom_line() +
  geom_point(size=2) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="Décès simulés", x="Temps (mois)", y="Nombre de décès", color="État") +
  theme_minimal()

p9 <- ggplot(donnees_sim, aes(x=temps, y=temperature, color=factor(etat))) +
  geom_line() +
  geom_point(size=2) +
  scale_color_manual(values=c("1"="orange", "2"="blue"), 
                     labels=c("1"="État 1", "2"="État 2")) +
  labs(title="Température simulée", x="Temps (mois)", y="Température (°C)", color="État") +
  theme_minimal()

grid.arrange(p7, p8, p9, ncol=1)

#------------------------------------------------
# PARTIE 4: COMPARAISON RÉEL VS SIMULÉ
#------------------------------------------------

# Création d'un dataframe combiné pour comparaison
donnees_reelles$type <- "Réel"
donnees_sim$type <- "Simulé"
donnees_combinees <- rbind(
  donnees_reelles[, c("temps", "deces", "temperature", "etat", "type")],
  transform(donnees_sim[, c("temps", "deces", "temperature", "etat", "type")], 
            temps = temps + max(donnees_reelles$temps))
)

# Visualisation comparative
p10 <- ggplot(donnees_combinees, aes(x=temps, y=deces, color=type)) +
  geom_line() +
  geom_point(size=2) +
  labs(title="Comparaison des décès: réel vs simulé", 
       x="Temps (mois)", y="Nombre de décès", color="Type de données") +
  theme_minimal() +
  geom_vline(xintercept=max(donnees_reelles$temps), linetype="dashed", color="grey50")

p11 <- ggplot(donnees_combinees, aes(x=temps, y=temperature, color=type)) +
  geom_line() +
  geom_point(size=2) +
  labs(title="Comparaison des températures: réel vs simulé", 
       x="Temps (mois)", y="Température (°C)", color="Type de données") +
  theme_minimal() +
  geom_vline(xintercept=max(donnees_reelles$temps), linetype="dashed", color="grey50")

grid.arrange(p10, p11, ncol=1)

# Distribution des valeurs réelles vs simulées
p12 <- ggplot(donnees_combinees, aes(x=deces, fill=type)) +
  geom_density(alpha=0.5) +
  labs(title="Distribution des décès", x="Nombre de décès", y="Densité", fill="Type") +
  theme_minimal()

p13 <- ggplot(donnees_combinees, aes(x=temperature, fill=type)) +
  geom_density(alpha=0.5) +
  labs(title="Distribution des températures", x="Température (°C)", y="Densité", fill="Type") +
  theme_minimal()

grid.arrange(p12, p13, ncol=1)

#------------------------------------------------
# PARTIE 5: PRÉDICTION DE SCÉNARIOS ALTERNATIFS
#------------------------------------------------

# Création d'un scénario avec réchauffement climatique
# (Augmentation progressive de la température)
donnees_sim_rechauffement <- donnees_sim
donnees_sim_rechauffement$temperature <- donnees_sim_rechauffement$temperature + 
  seq(0, 2, length.out=n)  # Augmentation progressive de 0 à 2°C

# Réajustement du modèle avec ces nouvelles données pour prédire l'impact sur les décès
modele_impact <- depmix(list(deces ~ 1, temperature ~ 1), 
                      data=donnees_sim_rechauffement,
                      nstates=nb_etats,
                      family=list(poisson(), gaussian()),
                      transition=~cos_t+sin_t)

# Utilisation des paramètres du modèle ajusté original
modele_impact <- setpars(modele_impact, getpars(modele_ajuste))

# Simulation avec le scénario de réchauffement
simulation_rechauffement <- simulate(modele_impact, nsim=1)

# Préparation pour visualisation
donnees_rechauffement <- data.frame(
  temps = 1:n,
  mois = rep(1:12, 3),
  annee = rep(2024:2026, each=12),
  deces = simulation_rechauffement$deces,
  temperature = donnees_sim_rechauffement$temperature,
  etat = simulation_rechauffement$state,
  scenario = "Réchauffement"
)

donnees_sim$scenario <- "Normal"
donnees_scenarios <- rbind(
  donnees_sim[, c("temps", "deces", "temperature", "etat", "scenario")],
  donnees_rechauffement[, c("temps", "deces", "temperature", "etat", "scenario")]
)

# Visualisation comparative des scénarios
p14 <- ggplot(donnees_scenarios, aes(x=temps, y=temperature, color=scenario)) +
  geom_line() +
  geom_point(size=2) +
  labs(title="Comparaison des températures selon les scénarios", 
       x="Temps (mois)", y="Température (°C)", color="Scénario") +
  theme_minimal()

p15 <- ggplot(donnees_scenarios, aes(x=temps, y=deces, color=scenario)) +
  geom_line() +
  geom_point(size=2) +
  labs(title="Impact du réchauffement sur les décès", 
       x="Temps (mois)", y="Nombre de décès", color="Scénario") +
  theme_minimal()

p16 <- ggplot(donnees_scenarios, aes(x=factor(etat), y=deces, fill=scenario)) +
  geom_boxplot() +
  labs(title="Distribution des décès par état et scénario", 
       x="État", y="Nombre de décès", fill="Scénario") +
  theme_minimal()

grid.arrange(p14, p15, p16, ncol=1)

# Résumé statistique des différences entre scénarios
diff_deces <- sum(donnees_rechauffement$deces) - sum(donnees_sim$deces)
pourcentage_diff <- (diff_deces / sum(donnees_sim$deces)) * 100

cat("Résumé de l'impact du réchauffement:\n")
cat("Différence totale de décès:", diff_deces, "\n")
cat("Pourcentage de changement:", round(pourcentage_diff, 2), "%\n")