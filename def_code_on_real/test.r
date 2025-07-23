# Chargement des packages nécessaires
library(depmixS4)
library(readxl)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)  # Pour gérer plusieurs couleurs d'états
library(grid)  # Pour utiliser grid.newpage()

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
trig_covs <- generate_trig_covariates(data$No_week, data$No_year, period, degree_trans_pol)

df = data.frame(
  death = data$Death_counts,
  death_rate = data$Death_counts / data$Weekly_exposure,
  temp = data$Temp_extrem,
  exposure = data$Weekly_exposure,
  log_exposure = log(data$Weekly_exposure),
  trend = data$No_year,
  trig_covs
)

df <- na.omit(df)
df <- df[is.finite(rowSums(df)), ]

# Ensure no NA, NaN, or Inf values in the covariates or response variables
df <- df[complete.cases(df), ]

#obs_formula_death_rate <- as.formula(paste("death_rate ~", paste(c(names(trig_covs)[-1], "offset(log_exposure)"), collapse = " + ")))
obs_formula_death_rate <- as.formula(paste("death_rate ~", paste(c(names(trig_covs)[-1]), collapse = " + ")))
obs_formula_temp <- as.formula(paste("temp ~", paste(c(names(trig_covs)[-1]), collapse = " + ")))

transition_formula <- as.formula(paste("~", paste(names(trig_covs)[c(-2, -1)], collapse = " + ")))

#------------------------------------------------
# FONCTION PRINCIPALE: ESTIMATION ET SIMULATION HMM
#------------------------------------------------

estimer_simuler_hmm <- function(donnees, nb_etats = 2, n_simulations = 1, 
                               seed = 123, maxit = 500, tol = 1e-10,
                               plot_results = TRUE) {

  
  
  # Définition du modèle HMM
  cat(paste("Définition d'un modèle HMM à", nb_etats, "états...\n"))
  modele_hmm <- depmix(list(obs_formula_death_rate, obs_formula_temp), 
                      data = df, #donnees,
                      nstates = nb_etats,
                      family = list(gaussian(link = log), gaussian()),
                      transition = transition_formula)
  
  # Affichage de la structure du modèle
  cat("Structure du modèle avant ajustement:\n")
  print(summary(modele_hmm))
  
  # Ajustement du modèle (estimation des paramètres)
  cat("Ajustement du modèle en cours, veuillez patienter...\n")
  set.seed(seed)  # Pour reproductibilité
  best_modele_ajuste <- NULL
  best_logLik <- -Inf

  for (i in 1:20) {  # Ajuster 20 fois
    modele_ajuste <- fit(modele_hmm, emcontrol = em.control(maxit = maxit, tol = tol))
    current_logLik <- logLik(modele_ajuste)
  
    if (current_logLik > best_logLik) {
      best_logLik <- current_logLik
      best_modele_ajuste <- modele_ajuste
    }
  }

  modele_ajuste <- best_modele_ajuste
  # Vérification de la convergence
  cat("Convergence du modèle:", ifelse(modele_ajuste@message == "Log likelihood converged to within tol. (relative change)", 
                                      "Réussie", paste("Échouée:", modele_ajuste@message)), "\n")
  
  # Affichage des paramètres estimés
  cat("Paramètres estimés:\n")
  print(summary(modele_ajuste))
  
  # Calcul des états les plus probables pour les données d'entrée
  etats_predits <- posterior(modele_ajuste, type = 'viterbi')
  donnees = df
  donnees$etat <- etats_predits$state
  
  # Utiliser le death_rate calibré pour calculer le death_count
  #donnees$death_rate_calibrated <- fitted(modele_ajuste)[, 1]
  #donnees$death_count_calculated <- rpois(nrow(donnees), lambda = donnees$Weekly_exposure * donnees$death_rate_calibrated)
  
  # Simulation de nouvelles données basées sur le modèle ajusté
  cat(paste("Simulation de", n_simulations, "jeux de données...\n"))
  set.seed(seed + 1)  # Différent du seed d'ajustement
  donnees_simulees_list <- list()
  
  for(i in 1:n_simulations) {
    sim_data <- simulate(modele_ajuste, nsim = 1)
    donnees$death_count_calculated 
    
    # Création d'un dataframe pour les données simulées
    n <- nrow(donnees)
    donnees_sim <- data.frame(
      time = 1:n,
      rate = sim_data@response[[1]][[1]]@y,
      death = rpois(nrow(donnees), lambda = donnees$exposure * sim_data@response[[1]][[1]]@y),
      temp = sim_data@response[[1]][[2]]@y,
      etat = donnees$etat,
      trig_covs,
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
    y_limits <- range(c(donnees$death, sim1$death))
    p9 <- ggplot(donnees, aes(x = factor(etat), y = death, fill = factor(etat))) +
      geom_boxplot() +
      scale_fill_manual(values = etat_palette) +
      labs(title = "Distribution des décès par état (données réelles)", 
           x = "État", y = "Nombre de décès", fill = "État") +
      theme_minimal() +
      ylim(y_limits)
    
    p10 <- ggplot(sim1, aes(x = factor(etat), y = death, fill = factor(etat))) +
      geom_boxplot() +
      scale_fill_manual(values = etat_palette) +
      labs(title = "Distribution des décès par état (données simulées)", 
           x = "État", y = "Nombre de décès", fill = "État") +
      theme_minimal() +
      ylim(y_limits)
    
    # Affichage des graphiques en grille
    cat("Affichage des visualisations...\n")
    png("visualization1.png", width=1500, height=1000)
    pdf("visualization1.pdf", width=1500, height=1000)
    grid.arrange(p1, p2, p3, ncol = 1)
    dev.off()

    png("visualization2.png", width=1500, height=1000)
    pdf("visualization2.pdf", width=1500, height=1000)
    grid.arrange(p4, p5, p6, ncol = 1)
    dev.off()

    png("visualization3.png", width=1500, height=1000)
    pdf("visualization3.pdf", width=1500, height=1000)
    grid.arrange(p7, p8, ncol = 2)
    dev.off()
    png("visualization4.png", width=1500, height=1000)
    pdf("visualization4.pdf", width=1500, height=1000)
    grid.arrange(p9, p10, ncol = 2)
    dev.off()
    # grid.newpage()
    # grid.arrange(p9, p10, ncol = 2)
    
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

# Comparaison de modèles avec différents nombres d'états
comparer_modeles <- function(donnees, max_etats = 4) {
  resultats_comparaison <- data.frame(
    nb_etats = 2:max_etats,
    logLik = NA,
    AIC = NA,
    BIC = NA,
    degree_of_freedom = NA,
    LRT_p_value = NA
  )
  
  for(i in 2:max_etats) {
    cat(paste("\n\n===== MODÈLE À", i, "ÉTATS =====\n\n"))
    res <- estimer_simuler_hmm(donnees, nb_etats = i, plot_results = FALSE)
    
    resultats_comparaison$logLik[i-1] <- res$log_likelihood[1]
    resultats_comparaison$AIC[i-1] <- res$AIC
    resultats_comparaison$BIC[i-1] <- res$BIC
    resultats_comparaison$degree_of_freedom[i-1] <- attr(logLik(res$modele_ajuste), "df")
  }
  # Calcul du test du rapport de vraisemblance (LRT)
  for(i in 2:max_etats) {
    lrt_stat <- 2 * (resultats_comparaison$logLik[i] - resultats_comparaison$logLik[i-1])
    lrt_df <- resultats_comparaison$degree_of_freedom[i] - resultats_comparaison$degree_of_freedom[i-1]
    lrt_p_value <- 1 - pchisq(lrt_stat, lrt_df)
    resultats_comparaison$LRT_p_value[i-1] <- lrt_p_value
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
resultats_comparaison <- comparer_modeles(df, max_etats = 3)
# print(resultats_comparaison)

# 2. Modèle avec nombre d'états spécifique (par exemple, le meilleur selon BIC)
resultats <- estimer_simuler_hmm(df, nb_etats = 3, n_simulations = 1, plot_results = TRUE)

# 3. Analyse du modèle optimal
# print(summary(resultats$modele_ajuste))

# COMMENTAIRE: Décommentez les lignes ci-dessus selon vos besoins
# Pour l'exécution complète, exécutez les 3 étapes
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

p3 <- ggplot(df, aes(x=time, y=exposure)) +
  geom_line() +
  geom_point() +
  labs(title="Exposure", x="Temps (semaine)", y="Exposure") +
  theme_minimal()


p4 <- ggplot(df, aes(x=temp, y=death)) +
  geom_point() +
  geom_smooth(method="loess") +
  labs(title="Relation entre température et décès", 
       x="Température (°C)", y="Nombre de décès") +
  theme_minimal()

grid.arrange(p1, p2, p3,p4, ncol=1)