# *****************************************************************************
#
#
# Purpose: Calibration of the 3-State Cancer Relative Survival (CRS) 
#          Markov Model using Nelder-Mead directed search algorithm
#
# Authors: 
# This work is developed by the Decision Analysis in R for Technologies in Health 
# (DARTH) workgroup:
#
# - Fernando Alarid-Escudero, PhD
# - Eva A. Enns, MS, PhD 
# - M.G. Myriam Hunink, MD, PhD 
# - Hawre J. Jalal, MD, PhD 
# - Eline Krijkamp, PhD 
# - Petros Pechlivanoglou, PhD
# - Alan Yang, MSc
#
# *****************************************************************************
#
# Notes:
#
# Please acknowledge our work. See details to cite below:
#
# - Alarid-Escudero F, MacLehose RF, Peralta Y, Kuntz KM, Enns EA. 
#   Non-identifiability in model calibration and implications for 
#   medical decision making. Med Decis Making. 2018; 38(7):810-821.
#
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, 
#   Hunink MG. An Overview of R in Health Decision Sciences. 
#   Med Decis Making. 2017; 37(3): 735-746.
#
# A walkthrough of the code could be found in the following link:
# - https://darth-git.github.io/calibSMDM2018-materials/
#
# *****************************************************************************

# ******************************************************************************
# 01 Calibration Overview ------------------------------------------------------
# ******************************************************************************

### 01.01 Model description  ---------------------------------------------------
# Model: 3-State Cancer Relative Survival (CRS) Markov Model
# Inputs to be calibrated: p_Mets, p_DieMets
# Targets: Survival data

### 01.02 Calibration method  --------------------------------------------------
# Search method: Directed search using Nelder-Mead algorithm
# Goodness-of-fit measure: Sum of Log-Likelihood

# ******************************************************************************
# 02 Setup ---------------------------------------------------------------------
# ******************************************************************************

### 02.01 Clear environment  ---------------------------------------------------
rm(list = ls())

### 02.02 Load packages  -------------------------------------------------------
# Install pacman if not present
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# Load pacman
library(pacman)

# Load (install if needed) CRAN packages
p_load(
  lhs,          # Latin Hypercube Sampling
  plotrix,      # Plotting with confidence intervals
  psych,        # Pairs panels
  mvtnorm,      # Multivariate normal distribution
  GGally,       # Pairwise plots
  ggplot2,      # Advanced plotting
  matrixcalc,   # To evaluate if matrix is positive definite
  Matrix        # To compute nearest positive definite matrix
)
# ******************************************************************************
# 03 Load calibration targets --------------------------------------------------
# ******************************************************************************

### 03.01 Load target data  ----------------------------------------------------
load("data/CRS_CalibTargets.RData")
lst_targets <- CRS_targets

### 03.02 Visualize calibration targets  ---------------------------------------
# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = lst_targets$Surv$time, 
                y = lst_targets$Surv$value, 
                ui = lst_targets$Surv$ub,
                li = lst_targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", 
                ylab = "Pr Survive")

# ******************************************************************************
# 04 Load model as a function --------------------------------------------------
# ******************************************************************************

### 04.01 Source model function  -----------------------------------------------
# Function inputs: parameters to be estimated through calibration
# Function outputs: model predictions corresponding to target data
source("code/Functions/CRS_MarkovModel_Function.R") # creates run_crs_markov()

### 04.02 Test model function  -------------------------------------------------
v_params_test <- c(p_Mets = 0.10, p_DieMets = 0.05)
run_crs_markov(v_params_test) # Test: function works correctly

# ******************************************************************************
# 05 Calibration specifications ------------------------------------------------
# ******************************************************************************

### 05.01 Set random seed  -----------------------------------------------------
set.seed(072218) # For reproducible sequence of random numbers

### 05.02 Define calibration parameters  ---------------------------------------
# Number of initial starting points
n_init <- 100

# Names and number of parameters to calibrate
v_param_names <- c("p_Mets", "p_DieMets")
n_param       <- length(v_param_names)

# Search space bounds
lb <- c(p_Mets = 0.04, p_DieMets = 0.04) # lower bound
ub <- c(p_Mets = 0.16, p_DieMets = 0.16) # upper bound

### 05.03 Define calibration targets  ------------------------------------------
v_target_names <- c("Surv")
n_target       <- length(v_target_names)

# ******************************************************************************
# 06 Goodness-of-fit function --------------------------------------------------
# ******************************************************************************

### 06.01 Define goodness-of-fit function  -------------------------------------
# Function to pass to Nelder-Mead algorithm
f_gof <- function(v_params){
  
  # Run model for parameter set v_params
  model_res <- run_crs_markov(v_params)
  
  # Calculate goodness-of-fit of model outputs to targets
  v_GOF <- numeric(n_target)
  
  # TARGET 1: Survival ("Surv")
  # Log likelihood  
  v_GOF[1] <- sum(dnorm(x = lst_targets$Surv$value,
                        mean = model_res$Surv,
                        sd = lst_targets$Surv$se,
                        log = T))
  
  # TARGET 2: (if you had more...)
  # log likelihood
  # v_GOF[2] <- sum(dnorm(x = lst_targets$Target2$value,
  #                        mean = model_res$Target2,
  #                        sd = lst_targets$Target2$se,
  #                        log = T))
  
  # OVERALL goodness-of-fit
  # Can give different targets different weights
  v_weights <- rep(1, n_target)
  
  # Weighted sum
  GOF_overall <- sum(v_GOF[1:n_target] * v_weights)
  
  # Return GOF
  return(GOF_overall)
}

# ******************************************************************************
# 07 Run calibration using Nelder-Mead -----------------------------------------
# ******************************************************************************

### 07.01 Record start time  ---------------------------------------------------
t_init <- Sys.time()

### 07.02 Sample initial starting values  --------------------------------------
# Generate random starting points for Nelder-Mead
v_params_init <- matrix(nrow = n_init, ncol = n_param)

for (i in 1:n_param) {
  v_params_init[, i] <- runif(n_init, min = lb[i], max = ub[i])
}

colnames(v_params_init) <- v_param_names

### 07.03 Initialize calibration results storage  ------------------------------
m_calib_res <- matrix(nrow = n_init, 
                      ncol = n_param + 1, 
                      dimnames = list(paste0("par_id_", 1:n_init),
                                      c(v_param_names, "Overall_fit")))

l_fit_nm <- vector(mode = "list", length = n_init)
names(l_fit_nm) <- paste0("par_id_", 1:n_init)

### 07.04 Run Nelder-Mead for each starting point  -----------------------------
for (j in 1:n_init) { # j <- 1
  
  ### use optim() as Nelder-Mead ###
  l_fit_nm[[j]] <- optim(v_params_init[j, ], f_gof,
                         control = list(fnscale = -1, # switches from minimization to maximization
                                        maxit = 1000), hessian = T)
  m_calib_res[j, ] <- c(l_fit_nm[[j]]$par, l_fit_nm[[j]]$value)
  
  ### to use a simulated annealing instead ###
  # fit_sa <- optim(v_params_init[j,], f_gof,
  #                method = c("SANN"),  # switches to using simulated annealing
  #                control = list(temp = 10, tmax = 10, # algorithm tuning parameters
  #                               fnscale = -1, maxit = 1000),
  #                hessian = T)
  # m_calib_res[j,] = c(fit_sa$par,fit_sa$value)
  
  ### to use a genetic algorithm instead ###
  # library(DEoptim)
  # f_fitness <- function(params){
  #   names(params) = v_param_names
  #   return(-f_gof(params))}
  # fit_ga = DEoptim(f_fitness, lower=lb, upper=ub)
  # m_calib_res[j,] = c(fit_ga$optim$bestmem,-1*fit_ga$optim$bestval)
  
}

### 07.05 Calculate computation time  ------------------------------------------
comp_time <- Sys.time() - t_init
comp_time

# ******************************************************************************
# 08 Explore best-fitting parameter sets ---------------------------------------
# ******************************************************************************

### 08.01 Sort results by goodness-of-fit  -------------------------------------
# Arrange parameter sets in order of fit
m_calib_res <- m_calib_res[order(-m_calib_res[, "Overall_fit"]), ]

# Best set
v_param_best <- m_calib_res[1, -3]

# Obtain id for best set
id_best_set <- rownames(m_calib_res)[1]

### 08.02 Examine top-performing parameter sets  -------------------------------
# Examine the top 10 best-fitting sets
m_calib_res[1:10, ]

### 08.03 Visualize top-performing parameter sets  -----------------------------
# Plot the top 10 (top 10%)
plot(m_calib_res[1:10, 1], m_calib_res[1:10, 2],
     xlim = c(lb[1], ub[1]), ylim = c(lb[2], ub[2]),
     xlab = colnames(m_calib_res)[1], ylab = colnames(m_calib_res)[2])

# Pairwise comparison of top 10 sets
pairs.panels(m_calib_res[1:10, v_param_names])

### 08.04 Compare best-fit model output to targets  ----------------------------
### Plot model-predicted output at mean vs targets ###
v_out_best <- run_crs_markov(m_calib_res[1, ])

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = lst_targets$Surv$time, y = lst_targets$Surv$value, 
                ui = lst_targets$Surv$ub,
                li = lst_targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(x = lst_targets$Surv$time, 
       y = v_out_best$Surv, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

# TARGET 2: (if you had more...)
# plotrix::plotCI(x = lst_targets$Target2$time, y = lst_targets$Target2$value, 
#                 ui = lst_targets$Target2$ub,
#                 li = lst_targets$Target2$lb,
#                 ylim = c(0, 1), 
#                 xlab = "Time", ylab = "Target 2")
# points(x = lst_targets$Target2$time, 
#        y = v_out_best$Target2, 
#        pch = 8, col = "red")
# legend("topright", 
#        legend = c("Target", "Model-predicted output"),
#        col = c("black", "red"), pch = c(1, 8))

# ******************************************************************************
# 09 Uncertainty quantification ------------------------------------------------
# ******************************************************************************

### 09.01 Hessian matrix  ------------------------------------------------------
m_hess <- l_fit_nm[[id_best_set]]$hessian
m_hess
# check if hessian is negative definite
eigen(m_hess)$values

### 09.02 Covariance matrix  ---------------------------------------------------
# Check if HESSIAN is Positive Definite; If not, make covariance 
# Positive Definite
# Is Positive Definite?
if (!is.positive.definite(-m_hess)) {
  print("Hessian is NOT Positive Definite")
  m_cov <- solve(-m_hess)
  print("Compute nearest positive definite matrix for COV matrix using `nearPD` function")
  m_cov <- Matrix::nearPD(m_cov)$mat
} else{
  print("Hessian IS Positive Definite")
  print("No additional adjustment to COV matrix")
  m_cov <- solve(-m_hess)
}

### 09.03 Correlation matrix  ---------------------------------------------------
m_cor <- cov2cor(m_cov)
m_cor

### 09.04 Standard errors  -----------------------------------------------------
m_se <- sqrt(diag(m_cov))
m_se

### 09.05 95% confidence interval  ---------------------------------------------
m_confint <- cbind(v_param_best - 1.96 * m_se,
                   v_param_best + 1.96 * m_se)
colnames(m_confint) <- c("LB", "UB")
m_confint

### 09.06 Draw sample of parameters  -------------------------------------------
n_samp <- 1000
m_param_best_sample <- rmvnorm(n = n_samp, 
                               mean = v_param_best, 
                               sigma = m_cov)
colnames(m_param_best_sample) <- v_param_names
pairs.panels(m_param_best_sample)

### 09.07 Fancier pairwise plot  -----------------------------------------------
gg_nm_pairs_corr <- GGally::ggpairs(data.frame(m_param_best_sample),
                                    upper = list(continuous = wrap("cor",
                                                                   color = "black",
                                                                   size = 5)),
                                    diag = list(continuous = wrap("barDiag",
                                                                  alpha = 0.8)),
                                    lower = list(continuous = wrap("points", 
                                                                   alpha = 0.3,
                                                                   size = 0.7)),
                                    columnLabels = v_param_names) +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "white",
                                        color = "white"),
        strip.text = element_text(hjust = 0))

gg_nm_pairs_corr
