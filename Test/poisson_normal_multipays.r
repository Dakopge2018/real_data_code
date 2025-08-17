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
path_to_data = "C:\\Users\\samue\\OneDrive\\Documents\\cours\\Projet de mémoire\\git\\real_data_code\\data\\main_database.xlsx"
OUTER_IMAGE = "C:\\Users\\samue\\OneDrive\\Documents\\cours\\Projet de mémoire\\git\\real_data_code\\image\\pn"

# Importation des données
data <- read_excel(path_to_data, sheet = "database_2")
#data <- read_excel("main_database.xlsx", sheet = "database")

# Check for and handle missing or infinite values and 0 values
data <- data[!is.na(data$Death_counts) & !is.infinite(data$Death_counts) & data$Death_counts > 0, ]
data <- na.omit(data)
data <- data[is.finite(rowSums(data[c('year','No_year', 'Temperature', 'Death_counts', 'Weekly_exposure', 'Temp_extrem')])), ]


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
  Code = data$Code,
  death = data$Death_counts,
  death_rate = data$Death_counts / data$Weekly_exposure,
  #temp = data$Temp_extrem,
  temp_norm = data$Temperature,
  temp = data$Temperature,
  #temp = data$Temp_new,
  exposure = data$Weekly_exposure,
  log_exposure = log(data$Weekly_exposure),
  trend = data$No_year,
  trig_covs
)

df <- na.omit(df)
# df <- df[is.finite(rowSums(df)), ]

df$date <- data$Date

# Ensure no NA, NaN, or Inf values in the covariates or response variables
df <- df[complete.cases(df), ]
donnees <- df
#obs_formula_death_rate <- as.formula(paste("death_rate ~", paste(c(names(trig_covs)[-1], "offset(log_exposure)"), collapse = " + ")))
obs_formula_death_rate <- as.formula(paste("death_rate ~", paste(c(names(trig_covs)[-1]), collapse = " + ")))
obs_formula_temp <- as.formula(paste("temp ~", paste(c(names(trig_covs)[-1]), collapse = " + ")))

transition_formula <- as.formula(paste("~", paste(names(trig_covs)[c(-2, -1)], collapse = " + ")))

#------------------------------------------------
# FONCTION PRINCIPALE: ESTIMATION ET SIMULATION HMM
#------------------------------------------------

estimer_simuler_hmm <- function(donnees, ext, nb_etats = 2, n_simulations = 1, 
                               seed = 1234, maxit = 500, tol = 1e-10,
                               plot_results = TRUE, n_ajustements = 20) {

  # Définition du modèle HMM
  cat(paste("Définition d'un modèle HMM à", nb_etats, "états...\n"))
  modele_hmm <- depmix(list(obs_formula_death_rate, obs_formula_temp), 
                      data = donnees, #donnees,
                      nstates = nb_etats,
                      family = list(gaussian(link ='log'), gaussian()),
                      transition = transition_formula)
  
  # Affichage de la structure du modèle
  cat("Structure du modèle avant ajustement:\n")
  print(summary(modele_hmm))
  
  # Ajustement du modèle avec sélection BIC/AIC
  cat(paste("Ajustement du modèle avec", n_ajustements, "tentatives...\n"))
  set.seed(seed)  # Pour reproductibilité
  
  # Initialisation des variables pour stocker les meilleurs modèles
  best_bic <- Inf
  best_aic <- Inf
  best_loglik <- -Inf
  best_modele_bic <- NULL
  best_modele_aic <- NULL
  best_modele_loglik <- NULL
  
  # Stocker tous les résultats pour analyse
  resultats_ajustement <- data.frame(
    iteration = integer(),
    loglik = numeric(),
    aic = numeric(),
    bic = numeric(),
    converged = logical()
  )

  for (i in 1:n_ajustements) {
    tryCatch({
      modele_ajuste <- fit(modele_hmm, emcontrol = em.control(maxit = maxit, tol = tol))
      
      # Calculer les critères
      current_loglik <- logLik(modele_ajuste)
      current_aic <- AIC(modele_ajuste)
      current_bic <- BIC(modele_ajuste)
      converged <- modele_ajuste@message == "Log likelihood converged to within tol. (relative change)"
      
      # Stocker les résultats
      # ...existing code...
    resultats_ajustement <- rbind(
        resultats_ajustement,
        data.frame(
            iteration = i,
            loglik = as.numeric(current_loglik),
            aic = current_aic,
            bic = current_bic,
            converged = converged,
            modele = I(list(modele_ajuste)) # Stocke le modèle dans une colonne liste
            )
        )
# ...existing code...
      
      # Mettre à jour le meilleur modèle selon BIC (critère principal)
      if (current_bic < best_bic) {
        best_bic <- current_bic
        best_modele_bic <- modele_ajuste
      }
      
      # Mettre à jour le meilleur modèle selon AIC
      if (current_aic < best_aic) {
        best_aic <- current_aic
        best_modele_aic <- modele_ajuste
      }
      
      # Mettre à jour le meilleur modèle selon log-vraisemblance
      if (current_loglik > best_loglik) {
        best_loglik <- current_loglik
        best_modele_loglik <- modele_ajuste
      }
      
    }, error = function(e) {
      cat("Erreur lors de l'ajustement", i, ":", e$message, "\n")
    })
    
    # Affichage du progrès
    if (i %% 10 == 0) {
      cat("Itération", i, "- Meilleur BIC:", round(best_bic, 2), "\n")
    }
  }

  # Sélection du modèle final (BIC comme critère principal)
  modele_ajuste <- best_modele_bic
  
  # Vérification que nous avons un modèle valide
  if (is.null(modele_ajuste)) {
    stop("Aucun ajustement de modèle n'a réussi. Vérifiez vos données et paramètres.")
  }
  
  # Affichage des résultats de sélection
  cat("\n=== RÉSULTATS DE LA SÉLECTION DE MODÈLE ===\n")
  cat("Nombre d'ajustements réussis:", nrow(resultats_ajustement), "/", n_ajustements, "\n")
  cat("Nombre d'ajustements convergés:", sum(resultats_ajustement$converged, na.rm = TRUE), "\n")
  cat("Meilleur BIC:", round(best_bic, 2), "\n")
  cat("Meilleur AIC:", round(best_aic, 2), "\n")
  cat("Meilleure log-vraisemblance:", round(best_loglik, 2), "\n")
  
  # Vérification de la cohérence entre BIC et AIC
  if (!is.null(best_modele_bic) && !is.null(best_modele_aic)) {
    bic_aic_coherent <- abs(BIC(best_modele_bic) - BIC(best_modele_aic)) < 1e-6
    if (bic_aic_coherent) {
      cat("✓ BIC et AIC sélectionnent le même modèle\n")
    } else {
      cat("⚠ BIC et AIC sélectionnent des modèles différents\n")
      cat("BIC du modèle sélectionné par AIC:", round(BIC(best_modele_aic), 2), "\n")
      cat("Différence BIC:", round(BIC(best_modele_aic) - best_bic, 2), "\n")
    }
  }
  
  # Statistiques sur les ajustements
  if (nrow(resultats_ajustement) > 0) {
    cat("\nStatistiques sur les ajustements:\n")
    cat("BIC - Min:", round(min(resultats_ajustement$bic, na.rm = TRUE), 2), 
        "Max:", round(max(resultats_ajustement$bic, na.rm = TRUE), 2), "\n")
    cat("AIC - Min:", round(min(resultats_ajustement$aic, na.rm = TRUE), 2), 
        "Max:", round(max(resultats_ajustement$aic, na.rm = TRUE), 2), "\n")
  }
  
  # Vérification de la convergence du modèle sélectionné
  cat("\nConvergence du modèle sélectionné:", 
      ifelse(modele_ajuste@message == "Log likelihood converged to within tol. (relative change)", 
             "Réussie", paste("Échouée:", modele_ajuste@message)), "\n")
  
  # Affichage des paramètres estimés
  cat("\nParamètres estimés du modèle sélectionné:\n")
  print(summary(modele_ajuste))
  
  # Calcul des états les plus probables pour les données d'entrée
  etats_predits <- posterior(modele_ajuste, type = 'viterbi')
  #donnees = df
  donnees$etat <- etats_predits$state
  
  # Simulation de nouvelles données basées sur le modèle ajusté
  cat(paste("Simulation de", n_simulations, "jeux de données...\n"))
  set.seed(seed + 1)  # Différent du seed d'ajustement
  donnees_simulees_list <- list()
  
  # for(i in 1:n_simulations) {
  #   n <- nrow(donnees)
  #   sim_data <- simulate(modele_ajuste, nsim = n)
  #  # donnees$death_count_calculated  
  #   # Création d'un dataframe pour les données simulées
    
  #   donnees_sim <- data.frame(
  #     time = 1:n,
  #     rate = sim_data@response[[1]][[1]]@y,
  #     death = (1 - exp(-sim_data@response[[1]][[1]]@y)) * donnees$exposure,
  #     temp = sim_data@response[[1]][[2]]@y,
  #     # etat = donnees$etat,
  #     etat_sim = sim_data@states,
  #     trig_covs,
  #     simulation_id = i
  #   )
    
  #   donnees_simulees_list[[i]] <- donnees_sim
  # }
  
    # ...existing code...
  for(i in 1:n_simulations) {
    n <- nrow(donnees)
    sim_data <- simulate(modele_ajuste, nsim = n)
    
    donnees_sim <- data.frame(
      time = 1:n,
      rate = sim_data@response[[1]][[1]]@y[1:n],
      death = (1 - exp(-sim_data@response[[1]][[1]]@y[1:n])) * donnees$exposure,
      temp = sim_data@response[[1]][[2]]@y[1:n],
      etat = donnees$etat,
      etat_sim = sim_data@states[1:n],
      trig_covs[1:n, ],
      simulation_id = i
    )
    
    donnees_simulees_list[[i]] <- donnees_sim
  }
  # ...existing code...
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
    
    # Graphique de convergence des critères (si package ggplot2 disponible)
    if (require(ggplot2, quietly = TRUE) && nrow(resultats_ajustement) > 1) {
      p_convergence <- ggplot(resultats_ajustement, aes(x = iteration)) +
        geom_line(aes(y = bic, color = "BIC")) +
        geom_line(aes(y = aic, color = "AIC")) +
        geom_point(aes(y = bic, color = "BIC"), alpha = 0.6) +
        geom_point(aes(y = aic, color = "AIC"), alpha = 0.6) +
        labs(title = "Évolution des critères d'information lors des ajustements",
             x = "Itération", y = "Valeur du critère", color = "Critère") +
        theme_minimal()
      
      #png(file.path(OUTER_IMAGE, "convergence_criteres.png"), width=1200, height=800)
      #print(p_convergence)
      #dev.off()
      
      pdf(file.path(OUTER_IMAGE, "convergence_criteres.pdf"), width=10, height=6)
      print(p_convergence)
      dev.off()
    }
    
    # Visualisation des données originales avec états prédits
    p1 <- ggplot(donnees, aes(x = date, y = factor(etat), color = factor(etat))) +
      geom_point(size = 3) +
      scale_color_manual(values = etat_palette) +
      labs(title = "États cachés prédits (données originales)", 
           x = "time", y = "État", color = "État") +
      theme_minimal()
    
    p2 <- ggplot(donnees, aes(x = date, y = death, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Décès par état (données originales)", 
           x = "Time", y = "Nombre de décès", color = "État") +
      theme_minimal()
    
    p3 <- ggplot(donnees, aes(x = date, y = temp_norm, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Température par état (données originales)", 
           x = "Time", y = "Température (°C)", color = "État") +
      theme_minimal()
      
    p3_bis <- ggplot(donnees, aes(x = date, y = temp, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Température extreme par état (données originales)", 
           x = "Time", y = "Température (°C)", color = "État") +
      theme_minimal()  

      
    
    # Visualisation des données simulées (première simulation)
    sim1 <- subset(donnees_simulees, simulation_id == 1)
    
    p4 <- ggplot(sim1, aes(x = time, y = factor(etat_sim), color = factor(etat_sim))) +
      geom_point(size = 3) +
      scale_color_manual(values = etat_palette) +
      labs(title = "États cachés simulés", 
           x = "time", y = "État", color = "État") +
      theme_minimal()
    
    p5 <- ggplot(sim1, aes(x = time, y = death, color = factor(etat_sim))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Décès simulés", 
           x = "time", y = "Nombre de décès", color = "État") +
      theme_minimal()
    
    p6 <- ggplot(sim1, aes(x = time, y = temp, color = factor(etat_sim))) +
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
    # png(file.path(OUTER_IMAGE, f"visualization1.png_{ext}"), width=1500, height=1000)
    # grid.arrange(p1, p2, p3,p3_bis, ncol = 1)
    # dev.off()
    pdf(file.path(OUTER_IMAGE, paste0("pn_visualization_real_data_",nb_etats, ext, ".pdf")), width=8, height=6)
    grid.arrange(p1, p2, p3, p3_bis, ncol = 1)
    dev.off()

    # png(file.path(OUTER_IMAGE, "visualization2.png"), width=1500, height=1000)
    # grid.arrange(p4, p5, p6, ncol = 1)
    # dev.off()
    pdf(file.path(OUTER_IMAGE, paste0("pn_visualization2_",nb_etats, ext, ".pdf")), width=8, height=6)
    grid.arrange(p4, p5, p6, ncol = 1)
    dev.off()

    # png(file.path(OUTER_IMAGE, "visualization3.png"), width=1500, height=1000)
    # grid.arrange(p7, p8, ncol = 2)
    # dev.off()
    pdf(file.path(OUTER_IMAGE, paste0("pn_visualization3_",nb_etats, ext, ".pdf")), width=10, height=8)
    grid.arrange(p7, p8, ncol = 2)
    dev.off()

    # png(file.path(OUTER_IMAGE, "visualization4.png"), width=1500, height=1000)
    # grid.arrange(p9, p10, ncol = 2)
    # dev.off()
    pdf(file.path(OUTER_IMAGE, paste0("pn_visualization4_",nb_etats, ext, ".pdf")), width=10, height=8)
    grid.arrange(p9, p10, ncol = 2)
    dev.off()
    
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
  
  # Retour des résultats enrichis
  resultats <- list(
    modele_original = modele_hmm,
    modele_ajuste = modele_ajuste,
    modeles_alternatifs = list(
      meilleur_aic = best_modele_aic,
      meilleure_loglik = best_modele_loglik
    ),
    donnees_avec_etats = donnees,
    donnees_simulees = donnees_simulees,
    criteres_selection = list(
      log_likelihood = best_loglik,
      BIC = best_bic,
      AIC = best_aic
    ),
    historique_ajustements = resultats_ajustement,
    coherence_criteres = if(!is.null(best_modele_bic) && !is.null(best_modele_aic)) {
      abs(BIC(best_modele_bic) - BIC(best_modele_aic)) < 1e-6
    } else { NA }
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
#resultats_comparaison <- comparer_modeles(df, max_etats = 2)
# print(resultats_comparaison)

# 2. Modèle avec nombre d'états spécifique (par exemple, le meilleur selon BIC)
# ...existing code...
resultats = list()
for (code in unique(df$Code)) {
  df_sub <- subset(df, Code == code)
  resultats[[as.character(code)]] <- estimer_simuler_hmm(df_sub, ext = code, nb_etats = 2, n_simulations = 1, plot_results = TRUE)
}


resultats <- estimer_simuler_hmm(df_sub, nb_etats = 2, n_simulations = 1, plot_results = FALSE)


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
pdf(file.path(OUTER_IMAGE, "distribution_initiale.pdf"), width=10, height=8)
grid.arrange(p1, p2, p3, p4, ncol = 1)
dev.off()