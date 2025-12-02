# *****************************************************************************
#
#
# Purpose: Bayesian calibration of the 3-State Cancer Relative Survival (CRS) 
#          Markov Model using Incremental Mixture Importance Sampling (IMIS)
#          with parallel computing
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
# Method: Incremental Mixture Importance Sampling (IMIS)
# Goodness-of-fit measure: Sum of Log-Likelihood
# Computing: Parallel processing enabled

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
  IMIS,         # Incremental Mixture Importance Sampling
  matrixStats,  # Summary statistics
  plotrix,      # Plotting with confidence intervals
  psych,        # Pairs panels
  ggplot2,      # Advanced plotting
  GGally,       # Pairwise plots
  dplyr,        # Data manipulation
  doParallel    # Parallel computing
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
# Number of posterior samples desired
n_resamp <- 1000

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
# 06 Prior distribution functions ----------------------------------------------
# ******************************************************************************

### 06.01 Function to sample from prior  ---------------------------------------
sample_prior <- function(n_samp){
  # Generate Latin Hypercube Sample
  m_lhs_unit   <- randomLHS(n = n_samp, k = n_param)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  
  # Transform to parameter space using uniform distribution
  for (i in 1:n_param) {
    m_param_samp[, i] <- qunif(m_lhs_unit[, i],
                               min = lb[i],
                               max = ub[i])
  }
  return(m_param_samp)
}

### 06.02 Visualize prior samples  ---------------------------------------------
pairs.panels(sample_prior(1000))

### 06.03 Function to calculate log-prior  -------------------------------------
calc_log_prior <- function(v_params){
  if (is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_samp <- nrow(v_params)
  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  
  # Calculate log-prior using uniform distribution
  for (i in 1:n_param) {
    lprior <- lprior + dunif(v_params[, i],
                             min = lb[i],
                             max = ub[i], 
                             log = TRUE)
  }
  return(lprior)
}

### 06.04 Function to calculate prior  -----------------------------------------
calc_prior <- function(v_params) { 
  exp(calc_log_prior(v_params)) 
}

### 06.05 Test prior functions  ------------------------------------------------
calc_log_prior(v_params = v_params_test)
calc_log_prior(v_params = sample_prior(10))
calc_prior(v_params = v_params_test)
calc_prior(v_params = sample_prior(10))

# ******************************************************************************
# 07 Likelihood functions ------------------------------------------------------
# ******************************************************************************

### 07.01 Function to calculate log-likelihood  --------------------------------
calc_log_lik <- function(v_params){
  if (is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_samp <- nrow(v_params)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
  llik_overall <- numeric(n_samp)
  
  for (j in 1:n_samp) {
    jj <- tryCatch({ 
      # Run model for parameter set v_params
      model_res <- run_crs_markov(v_params[j, ])
      
      # Calculate log-likelihood for TARGET 1: Survival
      v_llik[j, 1] <- sum(dnorm(x = lst_targets$Surv$value,
                                mean = model_res$Surv,
                                sd = lst_targets$Surv$se,
                                log = TRUE))
      
      # Calculate overall log-likelihood
      llik_overall[j] <- sum(v_llik[j, ])
    }, error = function(e) NA) 
    
    if (is.na(jj)) { llik_overall[j] <- -Inf }
  }
  return(llik_overall)
}

### 07.02 Function to calculate likelihood  ------------------------------------
calc_likelihood <- function(v_params){ 
  exp(calc_log_lik(v_params)) 
}

### 07.03 Test likelihood functions  -------------------------------------------
calc_log_lik(v_params = v_params_test)
calc_log_lik(v_params = sample_prior(10))
calc_likelihood(v_params = v_params_test)
calc_likelihood(v_params = sample_prior(10))

# ******************************************************************************
# 08 Posterior distribution functions ------------------------------------------
# ******************************************************************************

### 08.01 Function to calculate log-posterior  ---------------------------------
calc_log_post <- function(v_params) { 
  lpost <- calc_log_prior(v_params) + calc_log_lik(v_params)
  return(lpost) 
}

### 08.02 Function to calculate posterior  -------------------------------------
calc_post <- function(v_params) { 
  exp(calc_log_post(v_params)) 
}

### 08.03 Test posterior functions  --------------------------------------------
calc_log_post(v_params = v_params_test)
calc_log_post(v_params = sample_prior(10))
calc_post(v_params = v_params_test)
calc_post(v_params = sample_prior(10))


#' Parallel evaluation of log-likelihood function for a sets of parameters
#'
#' \code{log_lik_par} computes a log-likelihood value for one (or multiple) 
#' parameter set(s) using parallel computation.
#'
#' @param v_params Vector (or matrix) of model parameters.
#' @return 
#' A scalar (or vector) with log-likelihood values.
log_lik_par <- function(v_params,
                        ...) { 
  if (is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  n_samp <- nrow(v_params)
  
  ### Get OS
  os <- get_os()
  
  no_cores <- parallel::detectCores() - 1
  
  print(paste0("Parallelized Likelihood calculations on ", os, 
               " using ", no_cores, " cores"))
  
  n_time_init_likpar <- Sys.time()
  
  if (os == "macosx") {
    # Initialize cluster object
    cl <- parallel::makeForkCluster(no_cores) 
    doParallel::registerDoParallel(cl)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c) %dopar% {
      calc_log_lik(v_params[i, ]) # i = 1
    }
    n_time_end_likpar <- Sys.time()
  }
  if (os == "windows") {
    # Initialize cluster object
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
    opts <- list(attachExportEnv = TRUE)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c,
                              .export = ls(globalenv()),
                              .packages = c(),
                              .options.snow = opts) %dopar% {
                                calc_log_lik(v_params[i, ])
                              }
    n_time_end_likpar <- Sys.time()
  }
  if (os == "linux") {
    # Initialize cluster object
    cl <- parallel::makeCluster(no_cores)
    doMC::registerDoMC(cl)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c) %dopar% {
      calc_log_lik(v_params[i, ])
    }
    n_time_end_likpar <- Sys.time()
  }
  
  parallel::stopCluster(cl)
  n_time_total_likpar <- difftime(n_time_end_likpar, n_time_init_likpar, 
                                  units = "hours")
  print(paste0("Runtime: ", round(n_time_total_likpar, 2), " hrs."))
  #-# Try this: # PO
  rm(cl)        # PO
  gc()          # PO
  #-#           # PO
  return(v_llk)
}

#' Likelihood
#'
#' \code{likelihood} computes a likelihood value for one (or multiple) 
#' parameter set(s).
#'
#' @param v_params Vector (or matrix) of model parameters. 
#' @return 
#' A scalar (or vector) with likelihood values.
likelihood <- function(v_params){ 
  v_like <- exp(log_lik_par(v_params)) 
  return(v_like)
}

#' Get operating system
#' 
#' @return 
#' A string with the operating system.
#' @export
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "MacOSX"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
# ******************************************************************************
# 09 Run Bayesian calibration using IMIS ---------------------------------------
# ******************************************************************************

### 09.01 Record start time  ---------------------------------------------------
t_init <- Sys.time()

### 09.02 Define IMIS functions  -----------------------------------------------
prior          <- calc_prior
# likelihood     <- calc_likelihood
sample.prior   <- sample_prior

### 09.03 Run IMIS algorithm  --------------------------------------------------
fit_imis <- IMIS(
  B = 1000,         # Incremental sample size at each iteration
  B.re = n_resamp,  # Desired posterior sample size
  number_k = 10,    # Maximum number of iterations
  D = 0             # Number of samples to be deleted at each iteration
)

### 09.04 Extract posterior samples  -------------------------------------------
m_calib_res <- fit_imis$resample

### 09.05 Calculate posterior diagnostics  -------------------------------------
# Calculate log-likelihood and posterior probability for each sample
m_calib_res <- cbind(
  m_calib_res, 
  "Overall_fit" = calc_log_lik(m_calib_res[, v_param_names]),
  "Posterior_prob" = calc_post(m_calib_res[, v_param_names])
)

# Normalize posterior probabilities
m_calib_res[, "Posterior_prob"] <- m_calib_res[, "Posterior_prob"] / 
  sum(m_calib_res[, "Posterior_prob"])

### 09.06 Calculate computation time  ------------------------------------------
comp_time <- Sys.time() - t_init
comp_time

# ******************************************************************************
# 10 Explore posterior distribution --------------------------------------------
# ******************************************************************************

### 10.01 Plot posterior samples  ----------------------------------------------
# Color points by posterior probability
v_post_color <- scales::rescale(m_calib_res[, "Posterior_prob"])

plot(m_calib_res[, v_param_names],
     xlim = c(lb[1], ub[1]), 
     ylim = c(lb[2], ub[2]),
     xlab = v_param_names[1], 
     ylab = v_param_names[2],
     col = scales::alpha("black", v_post_color))

# Add centers of Gaussian components
points(fit_imis$center, col = "red", pch = 8)
legend("topright", 
       c("Draws from posterior", "Center of Gaussian components"),
       col = c("black", "red"), 
       pch = c(1, 8))

### 10.02 Pairwise plots with marginal histograms  -----------------------------
pairs.panels(m_calib_res[, v_param_names])

### 10.03 Calculate posterior summary statistics  ------------------------------
# Posterior mean
v_calib_post_mean <- colMeans(m_calib_res[, v_param_names])
v_calib_post_mean

# Posterior median and 95% credible interval
m_calib_res_95cr <- colQuantiles(m_calib_res[, v_param_names], 
                                 probs = c(0.025, 0.5, 0.975))
m_calib_res_95cr

# Maximum-a-posteriori (MAP) parameter set
v_calib_map <- m_calib_res[which.max(m_calib_res[, "Posterior_prob"]), ]

### 10.04 Compare model predictions to targets  --------------------------------
# Run model at MAP and posterior mean
v_out_best      <- run_crs_markov(v_calib_map[v_param_names])
v_out_post_mean <- run_crs_markov(v_calib_post_mean)

# Plot: Target vs MAP vs posterior mean
par(mar = c(5, 4, 4, 4)) 

plotrix::plotCI(x = lst_targets$Surv$time, 
                y = lst_targets$Surv$value, 
                ui = lst_targets$Surv$ub,
                li = lst_targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", 
                ylab = "Pr Survive")

points(x = lst_targets$Surv$time, 
       y = v_out_best$Surv, 
       pch = 8, col = "red")

points(x = lst_targets$Surv$time, 
       y = v_out_post_mean$Surv, 
       pch = 8, col = "blue")

legend("topright", 
       legend = c("Target", 
                  "Model-predicted output at MAP",
                  "Model-predicted output at posterior mean"),
       col = c("black", "red", "blue"), 
       pch = c(1, 8))

dev.off()

### 10.05 Advanced visualization: pairwise correlations  -----------------------
gg_post_pairs_corr <- GGally::ggpairs(
  data.frame(m_calib_res[, v_param_names]),
  upper = list(continuous = wrap("cor",
                                 color = "black",
                                 size = 5)),
  diag = list(continuous = wrap("barDiag",
                                alpha = 0.8)),
  lower = list(continuous = wrap("points", 
                                 alpha = 0.3,
                                 size = 0.7)),
  columnLabels = v_param_names
) +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "white",
                                        color = "white"),
        strip.text = element_text(hjust = 0))
gg_post_pairs_corr

### 10.06 Prior vs posterior comparison  ---------------------------------------
# Sample from prior
m_samp_prior <- sample.prior(n_resamp)

# Prepare data for plotting
df_samp_prior <- reshape2::melt(
  cbind(PDF = "Prior", as.data.frame(m_samp_prior)), 
  variable.name = "Parameter"
)

df_samp_post_imis <- reshape2::melt(
  cbind(PDF = "Posterior IMIS", 
        as.data.frame(m_calib_res[, v_param_names])),
  variable.name = "Parameter"
)

df_samp_prior_post <- dplyr::bind_rows(df_samp_prior, df_samp_post_imis)

# Plot prior vs posterior
gg_prior_post_imis <- ggplot(df_samp_prior_post, 
                             aes(x = value, y = ..density.., fill = PDF)) +
  facet_wrap(~Parameter, scales = "free", 
             ncol = 4,
             labeller = label_parsed) +
  scale_x_continuous(n.breaks = 6) +
  geom_density(alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
gg_prior_post_imis

# ******************************************************************************
# 11 Propagate parameter uncertainty -------------------------------------------
# ******************************************************************************

### 11.01 Initialize output matrix  --------------------------------------------
m_out_surv <- matrix(NA, 
                     nrow = n_resamp, 
                     ncol = length(lst_targets$Surv$value))

### 11.02 Run model for each posterior sample  ---------------------------------
for (i in 1:n_resamp) {
  model_res_temp <- run_crs_markov(m_calib_res[i, v_param_names])
  m_out_surv[i, ] <- model_res_temp$Surv
  
  # Progress indicator
  if (i/100 == round(i/100, 0)) { 
    cat('\r', paste(i/n_resamp * 100, "% done", sep = ""))
  }
}

### 11.03 Calculate posterior-predicted mean  ----------------------------------
v_out_surv_post_mean <- colMeans(m_out_surv)

### 11.04 Calculate posterior-predicted 95% credible interval  -----------------
m_out_surv_95cri <- colQuantiles(m_out_surv, probs = c(0.025, 0.975))

### 11.05 Prepare data for plotting  -------------------------------------------
df_out_post <- data.frame(
  Type = "Model output",
  bind_cols(Outcome = "Survival", 
            time = lst_targets[[1]]$time,
            value = v_out_surv_post_mean,
            lb = m_out_surv_95cri[, 1],
            ub = m_out_surv_95cri[, 2])
)

df_out_post$Outcome <- ordered(df_out_post$Outcome, levels = c("Survival"))

### 11.06 Plot targets vs model-predicted output  ------------------------------
df_targets <- data.frame(
  cbind(Type = "Target", 
        Outcome = "Survival", 
        lst_targets[[1]])
)

df_targets$Outcome <- ordered(df_targets$Outcome, levels = c("Survival"))

ggplot(df_targets, aes(x = time, y = value, ymin = lb, ymax = ub)) +
  geom_errorbar() +
  geom_line(data = df_out_post, 
            aes(x = time, y = value), 
            col = "blue") +
  geom_ribbon(data = df_out_post,
              aes(ymin = lb, ymax = ub), 
              alpha = 0.4, 
              fill = "blue") +
  facet_wrap(~ Outcome, scales = "free_x") +
  scale_x_continuous("Time", n.breaks = 8) +
  scale_y_continuous("Proportion", n.breaks = 8) +
  theme_bw(base_size = 16)

