# *****************************************************************************
#
# Script: CRS_Calibration_RandomSearch.R
#
# Purpose: Calibration of the 3-State Cancer Relative Survival (CRS) 
#          Markov Model using random search with Latin Hypercube Sampling
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
# Search method: Random search using Latin Hypercube Sampling
# Goodness-of-fit measure: Sum of log-likelihoods

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
  psych         # Pairs panels
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

# TARGET 2: (if you had more...)
# plotrix::plotCI(x = lst_targets$Target2$time, y = lst_targets$Target2$value, 
#                 ui = lst_targets$Target2$ub,
#                 li = lst_targets$Target2$lb,
#                 ylim = c(0, 1), 
#                 xlab = "Time", ylab = "Target 2")

# ******************************************************************************
# 04 Load model as a function --------------------------------------------------
# ******************************************************************************

### 04.01 Source model function  -----------------------------------------------
# Function inputs: parameters to be estimated through calibration
# Function outputs: model predictions corresponding to target data
source("code/CRS model/CRS_MarkovModel_Function.R") # creates run_crs_markov()

### 04.02 Test model function  -------------------------------------------------
v_params_test <- c(p_Mets = 0.10, p_DieMets = 0.05)
run_crs_markov(v_params_test) # Test: function works correctly

# ******************************************************************************
# 05 Calibration specifications ------------------------------------------------
# ******************************************************************************

### 05.01 Set random seed  -----------------------------------------------------
set.seed(072218) # For reproducible sequence of random numbers

### 05.02 Define calibration parameters  ---------------------------------------
# Number of random samples
n_samp <- 1000

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
# 06 Run calibration using Latin Hypercube Sampling ----------------------------
# ******************************************************************************

### 06.01 Record start time  ---------------------------------------------------
t_init <- Sys.time()

### 06.02 Generate random sample of parameter values  --------------------------
# Sample unit Latin Hypercube
m_lhs_unit <- randomLHS(n_samp, n_param)

# Rescale to min/max of each parameter
m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
for (i in 1:n_param) {
  m_param_samp[, i] <- qunif(m_lhs_unit[, i],
                             min = lb[i],
                             max = ub[i])
}
colnames(m_param_samp) <- v_param_names

# View resulting parameter set samples
pairs.panels(m_param_samp)

### 06.03 Run model for each set of parameter values  --------------------------
# Initialize goodness-of-fit matrix
m_GOF <- matrix(nrow = n_samp, ncol = n_target)
colnames(m_GOF) <- paste0(v_target_names, "_fit")

# Loop through sampled sets of input values
for (j in 1:n_samp) { # j <- 1
  
  # Run model for a given parameter set
  model_res <- run_crs_markov(v_params = m_param_samp[j, ])
  
  # Calculate goodness-of-fit of model outputs to targets
  
  # TARGET 1: Survival ("Surv")
  # Log likelihood  
  m_GOF[j, 1] <- sum(dnorm(x = lst_targets$Surv$value,
                           mean = model_res$Surv,
                           sd = lst_targets$Surv$se,
                           log = T))
  
  # weighted sum of squared errors (alternative to log likelihood)
  # w <- 1/(lst_targets$Surv$se^2)
  # m_GOF[j,1] <- -sum(w*(lst_targets$Surv$value - v_res)^2)
  
  # TARGET 2: (if you had more...)
  # log likelihood
  # m_GOF[j,2] <- sum(dnorm(x = lst_targets$Target2$value,
  #                        mean = model_res$Target2,
  #                        sd = lst_targets$Target2$se,
  #                        log = T))
  
} # End loop over sampled parameter sets

### 06.04 Combine fits to different targets into single GOF  -------------------
# Can give different targets different weights
v_weights <- matrix(1, nrow = n_target, ncol = 1)

# Matrix multiplication to calculate weighted sum of each GOF matrix row
v_GOF_overall <- c(m_GOF %*% v_weights)

# Store in GOF matrix with column name "Overall"
m_GOF <- cbind(m_GOF, Overall_fit = v_GOF_overall)

### 06.05 Calculate computation time  ------------------------------------------
comp_time <- Sys.time() - t_init

# ******************************************************************************
# 07 Explore best-fitting parameter sets ---------------------------------------
# ******************************************************************************

### 07.01 Sort results by goodness-of-fit  -------------------------------------
# Arrange parameter sets in order of fit
m_calib_res <- cbind(m_param_samp, m_GOF)
m_calib_res <- m_calib_res[order(-m_calib_res[, "Overall_fit"]), ]

### 07.02 Examine top-performing parameter sets  -------------------------------
# Examine the top 10 best-fitting sets
m_calib_res[1:10, ]

### 07.03 Visualize top-performing parameter sets  -----------------------------
# Plot the top 100 (top 10%)
plot(m_calib_res[1:100, 1], m_calib_res[1:100, 2],
     xlim = c(lb[1], ub[1]), ylim = c(lb[2], ub[2]),
     xlab = colnames(m_calib_res)[1],
     ylab = colnames(m_calib_res)[2])

# Pairwise comparison of top 100 sets
pairs.panels(m_calib_res[1:100, v_param_names])

### 07.04 Compare best and worst fit model outputs to targets  -----------------
# Plot model-predicted output at best and worst set vs targets
v_out_best  <- run_crs_markov(m_calib_res[1, ])
v_out_worst <- run_crs_markov(m_calib_res[999, ])

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = lst_targets$Surv$time, y = lst_targets$Surv$value, 
                ui = lst_targets$Surv$ub,
                li = lst_targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(x = lst_targets$Surv$time, 
       y = v_out_best$Surv, 
       pch = 8, col = "red")
points(x = lst_targets$Surv$time, 
       y = v_out_worst$Surv, 
       pch = 8, col = "blue")
legend("topright", 
       legend = c("Target", 
                  "Model-predicted output best set", 
                  "Model-predicted output worst set"),
       col = c("black", "red", "blue"), pch = c(1, 8, 8))

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

