# Libraries
library(depmixS4)
library(ggplot2)
library(readxl)
# Model parameters
set.seed(123)
states <- 4
degree_obs_pol <- 2
degree_trans_pol <- 2
period <- 52


# Read the Excel file
# Replace 'path_to_your_file.xlsx' with the actual path to your Excel file
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
    obs_poisson = round(data$Death_counts/100),
    obs_normal = data$Temperature,
    log_exposure = log(data$Weekly_exposure),
    trend = data$No_year,
    trig_covs)

df <- na.omit(df)
df <- df[is.finite(rowSums(df)), ]

# Ensure no NA, NaN, or Inf values in the covariates or response variables
df <- df[complete.cases(df), ]

# Transition formula
transition_formula <- as.formula(paste("~", paste(names(trig_covs)[c(-1,-2)], collapse = " + ")))

# Model
mod <- depmix(
  list(obs_poisson ~ trend + sin_1 + cos_1 + sin_2 + cos_2 + offset(log_exposure), obs_normal ~ trend + sin_1 + cos_1 + sin_2 + cos_2 ),
  data = df,
  nstates = states,
  family = list(poisson(link = "log"), gaussian()),
  transition = transition_formula,
)

#fitted_mod = fit(mod, verbose = TRUE)
set.seed(1)

fitted_model  <- multistart(mod,
  nstarts = 10,  # 10 initialisations
  initIters = 10,  # 10 itérations EM pour chaque initialisation
  emcontrol = em.control(
    maxit = 500,  # Max 500 itérations EM
    tol = 1e-08,  # Tolérance pour convergence
    crit = "relative",  # Critère de convergence
    random.start = TRUE,  # Randomisation des initialisations
    classification = "Hard"  # Classification 
    ))
summary(fitted_mod, which = 'transition')
summary(fitted_mod, which = 'response')

