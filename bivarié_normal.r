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
# Replace 'path_to_your_file.xlsx' with the actual path to your Excel file
data <- read_excel("main_database.xlsx", sheet = "database")

# Check for and handle missing or infinite values
data <- na.omit(data)
data <- data[is.finite(rowSums(data)), ]

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
    obs_poisson = data$Death_counts,
    obs_normal = data$Temperature,
    trig_covs)

df <- na.omit(df)
df <- df[is.finite(rowSums(df)), ]

# Transition formula
transition_formula <- as.formula(paste("~", paste(names(trig_covs)[c(-1,-2)], collapse = " + ")))

# Model
mod <- depmix(
  list(obs_poisson ~ sin_1 + cos_1 + trend, obs_normal ~ sin_1 + cos_1 + trend),
  data = df,
  nstates = 2,
  family = list(gaussian(), gaussian()),
  transition = transition_formula,
)

fitted_mod = fit(mod, verbose = TRUE)
set.seed(1)


#fitted_model  <- multistart(mod,
#  nstarts = 10,  # 10 initialisations
#  initIters = 10,  # 10 itérations EM pour chaque initialisation
#  emcontrol = em.control(
#    maxit = 500,  # Max 500 itérations EM
#    tol = 1e-08,  # Tolérance pour convergence
#    crit = "relative",  # Critère de convergence
#    random.start = TRUE,  # Randomisation des initialisations
#    classification = "Hard"  # Classification 
    ))
summary(fitted_mod, which = 'transition')
summary(fitted_mod, which = 'response')
