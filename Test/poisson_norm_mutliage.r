#Dataset with temperature (not extreme temperatures)
library(readxl)
# Importation des données
path_to_data = "C:\\Users\\samue\\OneDrive\\Documents\\cours\\Projet de mémoire\\git\\real_data_code\\data\\main_database.xlsx"
OUTER_IMAGE = "C:\\Users\\samue\\OneDrive\\Documents\\cours\\Projet de mémoire\\git\\real_data_code\\image\\pn\\multiple_ages"

# Importation des données

data <- read_excel(path_to_data, sheet = "database_3")

# Model parameters
set.seed(123)
states <- 2
degree_obs_pol <- 1
degree_trans_pol <- 1
period <- 52

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
  temp = data$Temperature,
  temp_extreme = data$Temp_extrem,
  exposure = data$Weekly_exposure,
  log_exposure = log(data$Weekly_exposure),
  trend = data$No_year,
  Age_factor = as.factor(data$Age_group),
  date = data$Date,
  trig_covs
)


# Load packages and color palette
library(ggplot2)
library(scico)
theme_set(theme_bw())
library(hmmTMB)
pal <- hmmTMB:::hmmTMB_cols
library(ggplot2)
library(gridExtra)
library(RColorBrewer)  # Pour gérer plusieurs couleurs d'états
library(grid)



sauvegarder_graphiques <- function(plots, repertoire_sortie = getwd()) {
  if (!dir.exists(repertoire_sortie)) {
    dir.create(repertoire_sortie, recursive = TRUE)
  }
  
  # Sauvegarder le graphique de convergence
  if (!is.null(plots$convergence)) {
    ggsave(file.path(repertoire_sortie, "convergence_criteres.png"), 
           plots$convergence, width = 12, height = 8, dpi = 300)
    ggsave(file.path(repertoire_sortie, "convergence_criteres.pdf"), 
           plots$convergence, width = 10, height = 6)
  }
  
  # Sauvegarder les graphiques des données originales
  if (!is.null(plots$donnees_originales)) {
    png(file.path(repertoire_sortie, "donnees_originales.png"), 
        width = 1500, height = 1000)
    do.call(grid.arrange, c(plots$donnees_originales, list(ncol = 1)))
    dev.off()
    
    pdf(file.path(repertoire_sortie, "donnees_originales.pdf"), 
        width = 8, height = 12)
    do.call(grid.arrange, c(plots$donnees_originales, list(ncol = 1)))
    dev.off()
  }
  
  # Sauvegarder les graphiques des données simulées
  if (!is.null(plots$donnees_simulees)) {
    png(file.path(repertoire_sortie, "donnees_simulees.png"), 
        width = 1500, height = 1000)
    do.call(grid.arrange, c(plots$donnees_simulees, list(ncol = 1)))
    dev.off()
    
    pdf(file.path(repertoire_sortie, "donnees_simulees.pdf"), 
        width = 8, height = 12)
    do.call(grid.arrange, c(plots$donnees_simulees, list(ncol = 1)))
    dev.off()
  }
  
  # Sauvegarder les graphiques de comparaison
  if (!is.null(plots$comparaison)) {
    png(file.path(repertoire_sortie, "comparaison_reel_simule.png"), 
        width = 1500, height = 800)
    do.call(grid.arrange, c(plots$comparaison, list(ncol = 2)))
    dev.off()
    
    pdf(file.path(repertoire_sortie, "comparaison_reel_simule.pdf"), 
        width = 12, height = 6)
    do.call(grid.arrange, c(plots$comparaison, list(ncol = 2)))
    dev.off()
  }
  
  # Sauvegarder les graphiques de distributions
  if (!is.null(plots$distributions)) {
    png(file.path(repertoire_sortie, "distributions_par_etat.png"), 
        width = 1500, height = 800)
    do.call(grid.arrange, c(plots$distributions, list(ncol = 2)))
    dev.off()
    
    pdf(file.path(repertoire_sortie, "distributions_par_etat.pdf"), 
        width = 12, height = 6)
    do.call(grid.arrange, c(plots$distributions, list(ncol = 2)))
    dev.off()
  }
  
  cat("Graphiques sauvegardés dans:", repertoire_sortie, "\n")
}


ajuster_hmm_optimal <- function(donnees, 
                                 nb_etats = 2, 
                                 n_simulations = 1000,
                                 seed = 123, 
                                 maxit = 1000, 
                                 tol = 1e-6,
                                 plot_results = TRUE, 
                                 n_ajustements = 10,
                                 sauvegarder_plots = FALSE,
                                 repertoire_sortie = getwd()) {
  
  # Charger les librairies nécessaires
  library(dplyr)
  library(ggplot2)
  
  # Définir le seed global
  set.seed(seed)
  
  # Fonction pour calculer les paramètres initiaux
  calculer_parametres_initiaux <- function(data, n_states) {
    # Calculer les statistiques de base
    mean_death <- mean(data$death, na.rm = TRUE)
    mean_temp <- mean(data$temp_extreme, na.rm = TRUE)
    
    # Diviser les données de température en groupes
    quantiles_temp <- quantile(data$temp_extreme, 
                              probs = seq(0, 1, length.out = n_states + 1), 
                              na.rm = TRUE)
    
    # Calculer les paramètres pour chaque état
    death_rates <- numeric(n_states)
    temp_means <- numeric(n_states)
    temp_sds <- numeric(n_states)
    
    for (i in 1:n_states) {
      if (i == 1) {
        mask <- data$temp_extreme <= quantiles_temp[i + 1]
      } else if (i == n_states) {
        mask <- data$temp_extreme > quantiles_temp[i]
      } else {
        mask <- data$temp_extreme > quantiles_temp[i] & 
                data$temp_extreme <= quantiles_temp[i + 1]
      }
      
      # Paramètres pour les taux de mortalité (avec variation aléatoire)
      base_death_rate <- mean(data$death[mask], na.rm = TRUE)
      if (is.na(base_death_rate) || base_death_rate <= 0) {
        base_death_rate <- mean_death
      }
      death_rates[i] <- max(base_death_rate * runif(1, 0.8, 1.2), 1)
      
      # Paramètres pour les températures
      temp_means[i] <- mean(data$temp_extreme[mask], na.rm = TRUE)
      temp_sds[i] <- sd(data$temp_extreme[mask], na.rm = TRUE)
      
      # Valeurs par défaut si NA
      if (is.na(temp_means[i])) temp_means[i] <- mean_temp
      if (is.na(temp_sds[i]) || temp_sds[i] <= 0) temp_sds[i] <- sd(data$temp_extreme, na.rm = TRUE)
    }
    
    return(list(
      death = list(rate = death_rates),
      temp_extreme = list(mean = temp_means, sd = temp_sds)
    ))
  }
  
  # Fonction pour ajuster un modèle HMM
  ajuster_modele_unique <- function(data, n_states, par_init, iteration) {
    tryCatch({
      # Formules
      f_transition <- ~ cos_1 + sin_1
      f_obs <- list(
        death = list(rate = ~Age_factor + trend + sin_1 + cos_1 + offset(log_exposure)),
        temp_extreme = list(mean = ~sin_1 + cos_1 + trend)
      )
      
      # Distributions (assumées définies globalement)
      # dists doit être défini dans l'environnement global
      
      # Créer les objets MarkovChain et Observation
      hid <- MarkovChain$new(data = data, n_states = n_states, formula = f_transition)
      obs <- Observation$new(data = data,
                            dists = dists,  # Assumé défini globalement
                            n_states = n_states,
                            par = par_init,
                            formulas = f_obs)
      
      # Créer et ajuster le HMM
      hmm <- HMM$new(obs = obs, hid = hid)
      hmm$fit(silent = TRUE, maxit = maxit, tol = tol)
      
      # Calculer les critères de sélection
      loglik <- logLik(hmm)[1]
      # aic <- -2 * loglik + 2 * length(hmm$coeff_fe())
      # bic <- -2 * loglik + log(nrow(data)) * length(hmm$coeff_fe())
      aic <- AIC(hmm)
      bic <- BIC(hmm)
      return(list(
        modele = hmm,
        loglik = loglik,
        aic = aic,
        bic = bic,
        iteration = iteration,
        convergence = TRUE
      ))
      
    }, error = function(e) {
      cat("Erreur dans l'ajustement", iteration, ":", e$message, "\n")
      return(list(
        modele = NULL,
        loglik = -Inf,
        aic = Inf,
        bic = Inf,
        iteration = iteration,
        convergence = FALSE,
        erreur = e$message
      ))
    })
  }
  
  # Fonction de simulation après ajustement
  simuler_donnees <- function(modele_optimal, n_sim, data_template) {
    if (is.null(modele_optimal)) {
      cat("Aucun modèle valide pour la simulation\n")
      return(NULL)
    }
    tryCatch({
      # Créer un template de données pour la simulation avec les bonnes variables
      template_sim <- data_template[1:n_sim, ]
      
      # Si on a moins de lignes que nécessaire, répéter les données
      if (nrow(data_template) < n_sim) {
        indices_repeat <- rep(seq_len(nrow(data_template)), length.out = n_sim)
        template_sim <- data_template[indices_repeat, ]
      }
    
      simulations <- modele_optimal$simulate(n = n_sim, data = template_sim)
      return(simulations)
    }, error = function(e) {
      cat("Erreur dans la simulation:", e$message, "\n")
      return(NULL)
    })
  }
  
  # AJUSTEMENT PRINCIPAL
  cat("Début de l'ajustement avec", n_ajustements, "tentatives...\n")
  
  resultats <- list()
  
  # Barre de progression
  pb <- txtProgressBar(min = 0, max = n_ajustements, style = 3)
  
  for (i in 1:n_ajustements) {
    # Générer des paramètres initiaux légèrement différents pour chaque tentative
    set.seed(seed + i)
    par_init <- calculer_parametres_initiaux(donnees, nb_etats)
    
    # Ajuster le modèle
    resultat <- ajuster_modele_unique(donnees, nb_etats, par_init, i)
    resultats[[i]] <- resultat
    
    # Mettre à jour la barre de progression
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Filtrer les modèles convergés
  modeles_valides <- resultats[sapply(resultats, function(x) x$convergence)]
  
  if (length(modeles_valides) == 0) {
    stop("Aucun modèle n'a convergé. Essayez d'augmenter maxit ou de modifier les paramètres.")
  }
  
  cat("\n", length(modeles_valides), "modèles ont convergé sur", n_ajustements, "tentatives.\n")
  
  # Sélectionner le meilleur modèle (plus faible AIC)
  aics <- sapply(modeles_valides, function(x) x$aic)
  indice_optimal <- which.min(aics)
  modele_optimal <- modeles_valides[[indice_optimal]]
  
  cat("Meilleur modèle: Itération", modele_optimal$iteration, "\n")
  cat("LogLik:", round(modele_optimal$loglik, 2), "\n")
  cat("AIC:", round(modele_optimal$aic, 2), "\n")
  cat("BIC:", round(modele_optimal$bic, 2), "\n")
  
  # Simulation
  cat("\nSimulation de", n_simulations, "observations...\n")
  df_2 <- donnees[, c("time", "sin_1", "cos_1", "trend", "exposure", "log_exposure","Age_factor")]
  donnees_simulees <- simuler_donnees(modele_optimal$modele, n_simulations, df_2)
  donnees_simulees$states <- attr(donnees_simulees, "state")  
  
 # Création des graphiques avancés
  plots <- list()
  
  if (plot_results && !is.null(donnees_simulees)) {
    # Charger les librairies nécessaires pour les visualisations
    if (!require(RColorBrewer, quietly = TRUE)) {
      install.packages("RColorBrewer")
      library(RColorBrewer)
    }
    if (!require(gridExtra, quietly = TRUE)) {
      install.packages("gridExtra")
      library(gridExtra)
    }
    
    # Créer un dossier pour les images (optionnel)
    if (!exists("OUTER_IMAGE")) {
      OUTER_IMAGE <- getwd()  # Utiliser le répertoire de travail actuel
    }
    
    # Obtenir les états prédits pour les données originales
    etats_predits <- modele_optimal$modele$viterbi()
    donnees$etat <- etats_predits
    
    # Créer une variable temporelle si elle n'existe pas
    if (!"date" %in% colnames(donnees) && !"time" %in% colnames(donnees)) {
      donnees$date <- seq_len(nrow(donnees))
    }
    if (!"time" %in% colnames(donnees)) {
      donnees$time <- donnees$date
    }
    
    # Préparer les données simulées
    if (is.data.frame(donnees_simulees)) {
      donnees_simulees$simulation_id <- 1
      donnees_simulees$states <- donnees_simulees$state  # Ajuster selon la structure
      donnees_simulees$time <- seq_len(nrow(donnees_simulees))
    }
    
    # Création d'une palette de couleurs adaptée au nombre d'états
    etat_palette <- brewer.pal(max(3, nb_etats), "Set1")
    if(nb_etats > 9) {
      # Interpolation si plus de 9 états
      etat_palette <- colorRampPalette(etat_palette)(nb_etats)
    }
    names(etat_palette) <- as.character(1:nb_etats)
    
    # Graphique de convergence des critères
    resultats_ajustement <- data.frame(
      iteration = sapply(modeles_valides, function(x) x$iteration),
      aic = sapply(modeles_valides, function(x) x$aic),
      bic = sapply(modeles_valides, function(x) x$bic),
      loglik = sapply(modeles_valides, function(x) x$loglik)
    )
    
    if (nrow(resultats_ajustement) > 1) {
      p_convergence <- ggplot(resultats_ajustement, aes(x = iteration)) +
        geom_line(aes(y = bic, color = "BIC")) +
        geom_line(aes(y = aic, color = "AIC")) +
        geom_point(aes(y = bic, color = "BIC"), alpha = 0.6) +
        geom_point(aes(y = aic, color = "AIC"), alpha = 0.6) +
        labs(title = "Évolution des critères d'information lors des ajustements",
             x = "Itération", y = "Valeur du critère", color = "Critère") +
        theme_minimal()
      
      plots$convergence <- p_convergence
      print(p_convergence)
    }
    
    # Visualisation des données originales avec états prédits
    p1 <- ggplot(donnees, aes(x = time, y = factor(etat), color = factor(etat))) +
      geom_point(size = 3) +
      scale_color_manual(values = etat_palette) +
      labs(title = "États cachés prédits (données originales)", 
           x = "Temps", y = "État", color = "État") +
      theme_minimal()
    
    p2 <- ggplot(donnees, aes(x = time, y = death, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Décès par état (données originales)", 
           x = "Temps", y = "Nombre de décès", color = "État") +
      theme_minimal()

    p3_bis <- ggplot(donnees, aes(x = time, y = temp, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Température par état (données originales)", 
           x = "Temps", y = "Température (°C)", color = "État") +
      theme_minimal()  
    
    p3 <- ggplot(donnees, aes(x = time, y = temp_extreme, color = factor(etat))) +
      geom_line() +
      geom_point(size = 2) +
      scale_color_manual(values = etat_palette) +
      labs(title = "Température extrême par état (données originales)", 
           x = "Temps", y = "Température (°C)", color = "État") +
      theme_minimal()
    
    plots$donnees_originales <- list(p1, p2, p3_bis, p3)
    
    # Visualisation des données simulées (si disponibles)
    if (is.data.frame(donnees_simulees)) {
      sim1 <- subset(donnees_simulees, simulation_id == 1)
      
      p4 <- ggplot(sim1, aes(x = time, y = factor(states), color = factor(states))) +
        geom_point(size = 3) +
        scale_color_manual(values = etat_palette) +
        labs(title = "États cachés simulés", 
             x = "Temps", y = "État", color = "État") +
        theme_minimal()
      
      p5 <- ggplot(sim1, aes(x = time, y = death, color = factor(states))) +
        geom_line() +
        geom_point(size = 2) +
        scale_color_manual(values = etat_palette) +
        labs(title = "Décès simulés", 
             x = "Temps", y = "Nombre de décès", color = "État") +
        theme_minimal()
      
      p6 <- ggplot(sim1, aes(x = time, y = temp_extreme, color = factor(states))) +
        geom_line() +
        geom_point(size = 2) +
        scale_color_manual(values = etat_palette) +
        labs(title = "Température simulée", 
             x = "Temps", y = "Température (°C)", color = "État") +
        theme_minimal()
      
      plots$donnees_simulees <- list(p4, p5, p6)
      
      # Comparaison entre réel et simulé
      donnees_temp <- donnees
      donnees_temp$type <- "Réel"
      sim1$type <- "Simulé"
      
      # S'assurer que les colonnes existent
      colonnes_communes <- intersect(c("time", "death", "temp_extreme", "type"), 
                                   intersect(colnames(donnees_temp), colnames(sim1)))
      
      if (length(colonnes_communes) >= 3) {
        donnees_comparaison <- rbind(
          donnees_temp[, colonnes_communes],
          sim1[, colonnes_communes]
        )
        
        p7 <- ggplot(donnees_comparaison, aes(x = time, y = death, color = type)) +
          geom_line() +
          geom_point(size = 2) +
          facet_wrap(~type) +
          labs(title = "Comparaison des décès: réel vs simulé", 
               x = "Temps", y = "Nombre de décès", color = "Type de données") +
          theme_minimal()
        
        p8 <- ggplot(donnees_comparaison, aes(x = time, y = temp_extreme, color = type)) +
          geom_line() +
          geom_point(size = 2) +
          facet_wrap(~type) +
          labs(title = "Comparaison des températures: réel vs simulé", 
               x = "Temps", y = "Température (°C)", color = "Type de données") +
          theme_minimal()
        
        plots$comparaison <- list(p7, p8)
      }
      
      # Distributions des valeurs par état
      if ("etat" %in% colnames(donnees) && "states" %in% colnames(sim1)) {
        y_limits <- range(c(donnees$death, sim1$death), na.rm = TRUE)
        
        p9 <- ggplot(donnees, aes(x = factor(etat), y = death, fill = factor(etat))) +
          geom_boxplot() +
          scale_fill_manual(values = etat_palette) +
          labs(title = "Distribution des décès par état (données réelles)", 
               x = "État", y = "Nombre de décès", fill = "État") +
          theme_minimal() +
          ylim(y_limits)
        
        p10 <- ggplot(sim1, aes(x = factor(states), y = death, fill = factor(states))) +
          geom_boxplot() +
          scale_fill_manual(values = etat_palette) +
          labs(title = "Distribution des décès par état (données simulées)", 
               x = "État", y = "Nombre de décès", fill = "État") +
          theme_minimal() +
          ylim(y_limits)
        
        plots$distributions <- list(p9, p10)
      }
    }
    
    # Affichage des graphiques principaux
    cat("Affichage des visualisations...\n")
    if (length(plots$donnees_originales) > 0) {
      print(do.call(grid.arrange, c(plots$donnees_originales, list(ncol = 1))))
    }
    
    if (length(plots$donnees_simulees) > 0) {
      print(do.call(grid.arrange, c(plots$donnees_simulees, list(ncol = 1))))
    }
    
    if (length(plots$comparaison) > 0) {
      print(do.call(grid.arrange, c(plots$comparaison, list(ncol = 2))))
    }
    
    if (length(plots$distributions) > 0) {
      print(do.call(grid.arrange, c(plots$distributions, list(ncol = 2))))
    }
    
    # Sauvegarder les graphiques si demandé
    if (sauvegarder_plots && length(plots) > 0) {
      sauvegarder_graphiques(plots, repertoire_sortie)
    }
  }
  
  # Retourner les résultats
  return(list(
    modele_optimal = modele_optimal$modele,
    tous_resultats = resultats,
    modeles_valides = modeles_valides,
    donnees_simulees = donnees_simulees,
    plots = plots,
    resume = list(
      n_ajustements = n_ajustements,
      n_converges = length(modeles_valides),
      meilleur_aic = modele_optimal$aic,
      meilleur_bic = modele_optimal$bic,
      meilleur_loglik = modele_optimal$loglik,
      iteration_optimale = modele_optimal$iteration
    )
  ))
}

# Exemple d'utilisation:
# Assurez-vous d'avoir défini 'dists' dans votre environnement global
dists <- list(death = "pois", temp_extreme = "norm")  # par exemple
# 
resultats <- ajuster_hmm_optimal(
  donnees = df,
  nb_etats = 2,
  n_simulations = 1000,
  seed = 123,
  maxit = 1000,
  tol = 1e-6,
  plot_results = TRUE,
  n_ajustements = 40,
  sauvegarder_plots = TRUE,  # Nouveau paramètre
  repertoire_sortie = OUTER_IMAGE
)

# Accéder au meilleur modèle
meilleur_modele <- resultats$modele_optimal

# Voir le résumé
print(resultats$resume)

# Coefficients du meilleur modèle
meilleur_modele$coeff_fe()
