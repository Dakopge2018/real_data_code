# Chargement des packages nécessaires
library(depmixS4)
library(readxl)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)  # Pour gérer plusieurs couleurs d'états

# Model parameters
set.seed(123)
states <- 4
degree_obs_pol <- 1
degree_trans_pol <- 1
period <- 52


#------------------------------------------------
# Importation des données

data <- read_excel("main_database.xlsx", sheet = "database")

# Check for and handle missing or infinite values
data <- na.omit(data)
data <- data[is.finite(rowSums(data[c('year','No_year', 'Temperature', 'Death_counts', 'Weekly_exposure')])), ]

generate_trig_covariates <- function(week_no, year_no, period, degree) {
  trig_covs <- data.frame(time = week_no, trend = year_no)
  for (d in 1:degree) {
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * week_no / period)
    trig_covs[[paste0("sin_", d)]] <- sin(2 * pi * d * week_no / period)
  }
  return(trig_covs)
}
trig_covs <- generate_trig_covariates(data$No_week, data$No_year,period, degree_trans_pol)

df = data.frame(
    death = round(data$Death_counts/100),
    temp = data$Temperature,
    log_exposure = log(data$Weekly_exposure),
    trend = data$No_year,
    trig_covs)

df <- na.omit(df)
df <- df[is.finite(rowSums(df)), ]

# Ensure no NA, NaN, or Inf values in the covariates or response variables
df <- df[complete.cases(df), ]
donnees <- df
#------------------------------------------------
# FONCTION PRINCIPALE: ESTIMATION ET SIMULATION HMM
#------------------------------------------------

estimer_simuler_hmm <- function(donnees, nb_etats = 2, n_simulations = 1, 
                               seed = 123, maxit = 500, tol = 1e-8,
                               plot_results = TRUE) {
  
  # Vérification des données d'entrée
  if(!all(c("death", "temp", "cos_1", "sin_1") %in% names(donnees))) {
    stop("Les données doivent contenir les colonnes: death, temp, cos_1, sin_1")
  }
  
  # Définition du modèle HMM
  cat(paste("Définition d'un modèle HMM à", nb_etats, "états...\n"))
  modele_hmm <- depmix(list(death~ sin_1 + cos_1 + trend + log_exposure, temp ~ sin_1 + cos_1 + trend), 
                      data = donnees,
                      nstates = nb_etats,
                      family = list(poisson(link = 'log'), gaussian()),
                      transition = ~ cos_1 + sin_1)
  
  # Affichage de la structure du modèle
  cat("Structure du modèle avant ajustement:\n")
  print(summary(modele_hmm))
  
  # Ajustement du modèle (estimation des paramètres)
  cat("Ajustement du modèle en cours, veuillez patienter...\n")
  set.seed(seed)  # Pour reproductibilité
  modele_ajuste <- fit(modele_hmm, emcontrol = em.control(maxit = maxit, tol = tol))
  
  # Vérification de la convergence
  cat("Convergence du modèle:", ifelse(modele_ajuste@message == "Log likelihood converged to within tol. (relative change)", 
                                      "Réussie", paste("Échouée:", modele_ajuste@message)), "\n")
  
  # Affichage des paramètres estimés
  cat("Paramètres estimés:\n")
  print(summary(modele_ajuste))
  
  # Calcul des états les plus probables pour les données d'entrée
  etats_predits <- posterior(modele_ajuste, type = 'viterbi')
  donnees$etat <- etats_predits$state
  
  # Simulation de nouvelles données basées sur le modèle ajusté
  cat(paste("Simulation de", n_simulations, "jeux de données...\n"))
  set.seed(seed + 1)  # Différent du seed d'ajustement
  donnees_simulees_list <- list()
  
  for(i in 1:n_simulations) {
    sim_data <- simulate(modele_ajuste, nsim = 1)
    
    # Création d'un dataframe pour les données simulées
    n <- nrow(donnees)
    donnees_sim <- data.frame(
      time = 1:n,
      death = sim_data@response[[1]][[1]]@y,
      temp = sim_data@response[[1]][[2]]@y,
      etat = donnees$etat,
      cos_1 = donnees$cos_1,  # Utilisation des mêmes covariables
      sin_1= donnees$sin_1,
      simulation_id = i
    )
    
    donnees_simulees_list[[i]] <- donnees_sim
  }
  
  # Combinaison des simulations en un seul dataframe
  donnees_simulees <- do.call(rbind, donnees_simulees_list)
  
  # Visualisation si demandée
  if(plot_results) {
    # Création d'une palette de couleurs adaptée au nombre d'états
    etat_palette <- brewer.pal(max(3, nb_etats), "Set1")
    if(nb_etats > 9) {
      # Interpolation si plus de 9 états
      etat_palette <- colorRampPalette(etat_palette)(nb_etats)
    }
    names(etat_palette) <- as.character(1:nb_etats)
    
    # Visualisation des données originales avec états prédits
    p1 <- ggplot(donnees, aes(x = time, y = factor(etat), color = factor(etat))) +
      geom_point(size = 3) +
      scale_color_manual(values = etat_palette) +
      labs(title = "États cachés prédits (données originales)", 
           x = "time", y = "État", color = "État") +
      theme_minimal()
    
    p2 <- ggplot(donnees, aes(x = time, y = death, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Décès par état (données originales)", 
           x = "Time", y = "Nombre de décès", color = "État") +
      theme_minimal()
    
    p3 <- ggplot(donnees, aes(x = time, y = temp, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Température par état (données originales)", 
           x = "Time", y = "Température (°C)", color = "État") +
      theme_minimal()
    
    # Visualisation des données simulées (première simulation)
    sim1 <- subset(donnees_simulees, simulation_id == 1)
    
    p4 <- ggplot(sim1, aes(x = time, y = factor(etat), color = factor(etat))) +
      geom_point(size = 3) +
      scale_color_manual(values = etat_palette) +
      labs(title = "États cachés simulés", 
           x = "time", y = "État", color = "État") +
      theme_minimal()
    
    p5 <- ggplot(sim1, aes(x = time, y = death, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Décès simulés", 
           x = "time", y = "Nombre de décès", color = "État") +
      theme_minimal()
    
    p6 <- ggplot(sim1, aes(x = time, y = temp, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Température simulée", 
           x = "time", y = "Température (°C)", color = "État") +
      theme_minimal()
    
    # Comparaison entre réel et simulé
    donnees$type <- "Réel"
    sim1$type <- "Simulé"
    donnees_comparaison <- rbind(
      donnees[, c("time", "death", "temp", "type")],
      sim1[, c("time", "death", "temp", "type")]
    )
    
    p7 <- ggplot(donnees_comparaison, aes(x = time, y = death, color = type)) +
      geom_line() +
      geom_point(size = 2) +
      facet_wrap(~type) +
      labs(title = "Comparaison des décès: réel vs simulé", 
           x = "time", y = "Nombre de décès", color = "Type de données") +
      theme_minimal()
    
    p8 <- ggplot(donnees_comparaison, aes(x = time, y = temp, color = type)) +
      geom_line() +
      geom_point(size = 2) +
      facet_wrap(~type) +
      labs(title = "Comparaison des températures: réel vs simulé", 
           x = "time", y = "Température (°C)", color = "Type de données") +
      theme_minimal()
    
    # Distributions des valeurs par état
    p9 <- ggplot(donnees, aes(x = factor(etat), y = death, fill = factor(etat))) +
      geom_boxplot() +
      scale_fill_manual(values = etat_palette) +
      labs(title = "Distribution des décès par état (données réelles)", 
           x = "État", y = "Nombre de décès", fill = "État") +
      theme_minimal()
    
    p10 <- ggplot(sim1, aes(x = factor(etat), y = death, fill = factor(etat))) +
      geom_boxplot() +
      scale_fill_manual(values = etat_palette) +
      labs(title = "Distribution des décès par état (données simulées)", 
           x = "État", y = "Nombre de décès", fill = "État") +
      theme_minimal()
    
    # Affichage des graphiques en grille
    cat("Affichage des visualisations...\n")
    grid.arrange(p1, p2, p3, ncol = 1)
    grid.arrange(p4, p5, p6, ncol = 1)
    grid.arrange(p7, p8, ncol = 2)
    grid.arrange(p9, p10, ncol = 2)
    
    # Si plusieurs simulations, visualiser la variabilité
    if(n_simulations > 1) {
      p_var <- ggplot(donnees_simulees, aes(x = time, y = death, group = interaction(simulation_id, etat), 
                                           color = factor(etat))) +
        geom_line(alpha = 0.3) +
        scale_color_manual(values = etat_palette) +
        labs(title = paste("Variabilité sur", n_simulations, "simulations"), 
             x = "time", y = "Nombre de décès", color = "État") +
        theme_minimal()
      
      print(p_var)
    }
  }
  
  # Retour des résultats
  resultats <- list(
    modele_original = modele_hmm,
    modele_ajuste = modele_ajuste,
    donnees_avec_etats = donnees,
    donnees_simulees = donnees_simulees,
    log_likelihood = logLik(modele_ajuste),
    BIC = BIC(modele_ajuste),
    AIC = AIC(modele_ajuste)
  )
  
  return(resultats)
}

#------------------------------------------------
# EXEMPLE D'UTILISATION
#------------------------------------------------

# Génération de données d'exemple (à remplacer par vos données réelles)
# generer_donnees_exemple <- function(n = 36, seed = 123) {
#   set.seed(seed)
#   time <- 1:n
  
#   # Covariables saisonnières
#   cos_t <- cos(2*pi*time/12)
#   sin_t <- sin(2*pi*time/12)
  
#   # Génération de variables avec différents états
#   # Imaginons 3 états: hiver rigoureux, saison intermédiaire, été chaud
#   etats_vrais <- rep(NA, n)
  
#   # Saisons: hiver (1), printemps/automne (2), été (3)
#   for(i in 1:n) {
#     mois <- ((i-1) %% 12) + 1
#     if(mois %in% c(12, 1, 2)) etats_vrais[i] <- 1  # Hiver
#     else if(mois %in% c(6, 7, 8)) etats_vrais[i] <- 3  # Été
#     else etats_vrais[i] <- 2  # Printemps/Automne
#   }
  
#   # Températures selon les états
#   temp_base <- c(3, 15, 25)  # Moyennes par état
#   temp_sd <- c(3, 4, 2)      # Écarts-types par état
  
#   temperature <- numeric(n)
#   for(i in 1:n) {
#     etat <- etats_vrais[i]
#     temperature[i] <- rnorm(1, temp_base[etat], temp_sd[etat])
#   }
  
#   # Décès selon les états
#   deces_lambda <- c(120, 90, 75)  # Moyennes par état
  
#   deces <- numeric(n)
#   for(i in 1:n) {
#     etat <- etats_vrais[i]
#     # Une légère dépendance négative à la température
#     lambda_ajuste <- deces_lambda[etat] * exp(-0.01 * (temperature[i] - temp_base[etat]))
#     deces[i] <- rpois(1, lambda_ajuste)
#   }
  
#   # Création du dataframe
#   data.frame(
#     time = time,
#     mois = rep(1:12, length.out = n),
#     annee = rep(2021:(2021+ceiling(n/12)-1), each = 12)[1:n],
#     cos_t = cos_t,
#     sin_t = sin_t,
#     deces = deces,
#     temperature = temperature,
#     etat_vrai = etats_vrais
#   )
# }

# Génération de données d'exemple
# donnees_exemple <- generer_donnees_exemple(n = 48)  # 4 ans de données

# Visualisation des données d'exemple
p1 <- ggplot(df, aes(x=time, y=death)) +
  geom_line() +
  geom_point() +
  labs(title="Décès observés", x="Temps (semaine)", y="Nombre de décès") +
  theme_minimal()

p2 <- ggplot(df, aes(x=time, y=temp)) +
  geom_line() +
  geom_point() +
  labs(title="Température observée", x="Temps (semaine)", y="Température (°C)") +
  theme_minimal()

p3 <- ggplot(df, aes(x=temp, y=death)) +
  geom_point() +
  geom_smooth(method="loess") +
  labs(title="Relation entre température et décès", 
       x="Température (°C)", y="Nombre de décès") +
  theme_minimal()

grid.arrange(p1, p2, p3, ncol=1)

# Comparaison de modèles avec différents nombres d'états
comparer_modeles <- function(donnees, max_etats = 4) {
  resultats_comparaison <- data.frame(
    nb_etats = 2:max_etats,
    logLik = NA,
    AIC = NA,
    BIC = NA
  )
  
  for(i in 2:max_etats) {
    cat(paste("\n\n===== MODÈLE À", i, "ÉTATS =====\n\n"))
    res <- estimer_simuler_hmm(donnees, nb_etats = i, plot_results = FALSE)
    
    resultats_comparaison$logLik[i] <- res$log_likelihood
    resultats_comparaison$AIC[i] <- res$AIC
    resultats_comparaison$BIC[i] <- res$BIC
  }
  
  # Graphiques de comparaison
  p1 <- ggplot(resultats_comparaison, aes(x = nb_etats, y = AIC)) +
    geom_line() +
    geom_point(size = 3) +
    labs(title = "Comparaison des modèles: AIC", 
         x = "Nombre d'états", y = "AIC (plus petit = meilleur)") +
    theme_minimal()
  
  p2 <- ggplot(resultats_comparaison, aes(x = nb_etats, y = BIC)) +
    geom_line() +
    geom_point(size = 3) +
    labs(title = "Comparaison des modèles: BIC", 
         x = "Nombre d'états", y = "BIC (plus petit = meilleur)") +
    theme_minimal()
  
  grid.arrange(p1, p2, ncol = 2)
  
  return(resultats_comparaison)
}

# Exemple d'utilisation:

# 1. Comparer différents nombres d'états
resultats_comparaison <- comparer_modeles(df, max_etats = 4)
print(resultats_comparaison)

# 2. Modèle avec nombre d'états spécifique (par exemple, le meilleur selon BIC)
resultats <- estimer_simuler_hmm(df, nb_etats = 6, n_simulations = 1)

# 3. Analyse du modèle optimal
# print(summary(resultats$modele_ajuste))

# COMMENTAIRE: Décommentez les lignes ci-dessus selon vos besoins
# Pour l'exécution complète, exécutez les 3 étapes