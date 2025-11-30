# *****************************************************************************
#
#
# Purpose: Deterministic and Probabilistic Sensitivity Analysis for the 
#          Sick-Sicker model
#
# Authors: 
# This work is developed by the Decision Analysis in R for Technologies in Health 
# (DARTH) workgroup:
#
# - Fernando Alarid-Escudero, PhD <falarid@stanford.edu>
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
# Please cite the article when using this code.
#
# This code implements a simulation-time-dependent Sick-Sicker cSTM model to 
# conduct a CEA of two strategies:
# - Standard of Care (SoC): best available care for the patients with the disease. 
#   This scenario reflects the natural history of the disease progression.
# - Strategy AB: This strategy combines treatment A and treatment B. The disease 
#   progression is reduced, and individuals in the Sick state have an improved 
#   quality of life.
#
# *****************************************************************************

# ******************************************************************************
# 01 Setup ---------------------------------------------------------------------
# ******************************************************************************

### 01.01 Clear environment  ---------------------------------------------------
rm(list = ls())      # clear memory (removes all the variables from the workspace)

### 01.02 Load packages  -------------------------------------------------------
# Install pacman if not present
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# Load pacman
library(pacman)

# Load (install if needed) CRAN packages
p_load(
  darthtools,  # DARTH tools
  dampack,     # Decision-analytic modeling package
  reshape2,    # Data reshaping
  dplyr,       # Data manipulation
  ggplot2,     # Advanced plotting
  ggthemes,    # Plot themes
  scales       # Scale functions for visualization
)

### 01.03 Load model functions  ------------------------------------------------
# Load model, CEA and PSA functions 
source('code/Functions_markov_sick-sicker_time.R')

# ******************************************************************************
# 02 Model input parameters ----------------------------------------------------
# ******************************************************************************

### 02.01 General setup  -------------------------------------------------------
cycle_length <- 1   # cycle length equal to one year (use 1/12 for monthly)
n_age_init   <- 25  # age at baseline
n_age_max    <- 100 # maximum age of follow up
n_cycles     <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles

# Age labels 
v_age_names  <- paste(rep(n_age_init:(n_age_max - 1), each = 1/cycle_length), 
                      1:(1/cycle_length), 
                      sep = ".")

### 02.02 Age-dependent mortality rates  ---------------------------------------
lt_usa_2015 <- read.csv("data/HMD_USA_Mx_2015.csv")

# Extract age-specific all-cause mortality for ages in model time horizon
v_r_mort_by_age <- lt_usa_2015 %>% 
  dplyr::filter(Age >= n_age_init & Age < n_age_max) %>%
  dplyr::select(Total) %>%
  as.matrix()

# Process model inputs 
# Age-specific transition rates to the Dead state for all cycles 
v_r_HDage  <- rep(v_r_mort_by_age, each = 1/cycle_length)

# Name age-specific mortality vector 
names(v_r_HDage) <- v_age_names

### 02.03 Strategies  ----------------------------------------------------------
v_names_str <- c("Standard of care", # store the strategy names
                 "Strategy AB") 
n_str       <- length(v_names_str)   # number of strategies

### 02.04 List of input parameters  --------------------------------------------
l_params_all <- list(
  # Transition probabilities (per cycle), hazard ratios
  v_r_HDage = v_r_HDage, # constant rate of dying when Healthy (all-cause mortality)
  r_HS1     = 0.15,      # constant annual rate of becoming Sick when Healthy conditional on surviving
  r_S1H     = 0.5,       # constant annual rate of becoming Healthy when Sick conditional on surviving
  r_S1S2    = 0.105,     # constant annual rate of becoming Sicker when Sick conditional on surviving
  hr_S1     = 3,         # hazard ratio of death in Sick vs Healthy 
  hr_S2     = 10,        # hazard ratio of death in Sicker vs Healthy 
  # Effectiveness of treatment AB 
  hr_S1S2_trtAB = 0.6,   # hazard ratio of becoming Sicker when Sick under treatment AB
  ## State rewards
  # Costs
  c_H       = 2000,      # cost of remaining one cycle in Healthy 
  c_S1      = 4000,      # cost of remaining one cycle in Sick 
  c_S2      = 15000,     # cost of remaining one cycle in Sicker 
  c_D       = 0,         # cost of being dead (per cycle)
  c_trtAB   = 25000,     # cost of treatment A
  # Utilities
  u_H       = 1,         # utility when Healthy 
  u_S1      = 0.75,      # utility when Sick 
  u_S2      = 0.5,       # utility when Sicker
  u_D       = 0,         # utility when Dead 
  u_trtAB   = 0.95,      # utility when being treated with A
  ## Transition rewards
  du_HS1    = 0.01,      # disutility when transitioning from Healthy to Sick
  ic_HS1    = 1000,      # increase in cost when transitioning from Healthy to Sick
  ic_D      = 2000,      # increase in cost when dying
  # Initial and maximum ages
  n_age_init = 25,
  n_age_max  = 100,
  # Discount rates
  d_c        = 0.03,     # annual discount rate for costs 
  d_e        = 0.03,     # annual discount rate for QALYs,
  # Cycle length
  cycle_length = 1
)

### 02.05 Test model calculation  ----------------------------------------------
calculate_ce_out(l_params_all = l_params_all, verbose = TRUE)

# ******************************************************************************
# 03 One-way sensitivity analysis (OWSA) ---------------------------------------
# ******************************************************************************

### 03.01 Define OWSA parameters  ----------------------------------------------
options(scipen = 999) # disabling scientific notation in R

# data.frame containing all parameters, their base-case values, and the min and 
# max values of the parameters of interest 
df_params_owsa <- data.frame(
  pars = c("r_S1S2", "c_trtAB", "u_S1", "u_trtAB"),
  min  = c(0.05 , 20000 , 0.65, 0.7),  # min parameter values
  max  = c(0.155, 40000 , 0.85, 0.95)  # max parameter values
)

### 03.02 Run OWSA  ------------------------------------------------------------
owsa_nmb  <- run_owsa_det(
  params_range    = df_params_owsa,   # data.frame with parameters for OWSA
  params_basecase = l_params_all,     # list with all parameters
  nsamp           = 100,              # number of parameter values
  FUN             = calculate_ce_out, # function to compute outputs
  outcomes        = c("NMB"),         # output to do the OWSA on
  strategies      = v_names_str,      # names of the strategies
  n_wtp           = 120000            # extra argument to pass to FUN
)

### 03.03 Plot OWSA results  ---------------------------------------------------
plot(owsa_nmb, txtsize = 10, n_x_ticks = 4, 
     facet_scales = "free") +
  theme(legend.position = "bottom")

### 03.04 Optimal strategy with OWSA  ------------------------------------------
owsa_opt_strat(owsa = owsa_nmb, txtsize = 10)

### 03.05 Tornado plot  --------------------------------------------------------
owsa_tornado(owsa = owsa_nmb)

# ******************************************************************************
# 04 Two-way sensitivity analysis (TWSA) ---------------------------------------
# ******************************************************************************

### 04.01 Define TWSA parameters  ----------------------------------------------
# dataframe containing all parameters, their basecase values, and the min and 
# max values of the parameters of interest
df_params_twsa <- data.frame(
  pars = c("c_trtAB", "u_trtAB"),
  min  = c(18000, 0.80),  # min parameter values
  max  = c(36000, 0.98)   # max parameter values
)

### 04.02 Run TWSA  ------------------------------------------------------------
twsa_nmb <- run_twsa_det(
  params_range    = df_params_twsa,   # data.frame with parameters for TWSA
  params_basecase = l_params_all,     # list with all parameters
  nsamp           = 40,               # number of parameter values
  FUN             = calculate_ce_out, # function to compute outputs
  outcomes        = c("NMB"),         # output to do the TWSA on
  strategies      = v_names_str,      # names of the strategies
  n_wtp           = 120000            # extra argument to pass to FUN
)

### 04.03 Plot TWSA  -----------------------------------------------------------
plot(twsa_nmb)

# ******************************************************************************
# 05 Probabilistic Sensitivity Analysis (PSA) ----------------------------------
# ******************************************************************************

### 05.01 Model input for PSA  -------------------------------------------------
# Store the parameter names into a vector
v_names_params <- names(l_params_all)

### 05.02 Test PSA functions  --------------------------------------------------
# Test function to compute CE outcomes
calculate_ce_out(l_params_all) 

# Test function to generate PSA input dataset
generate_psa_params(10) 

### 05.03 Generate PSA dataset  ------------------------------------------------
# Number of simulations
n_sim <- 1000

# Generate PSA input dataset
df_psa_input <- generate_psa_params(n_sim = n_sim)

# First six observations
head(df_psa_input)

### 05.04 Visualize parameter distributions  -----------------------------------
# Histogram of parameters 
ggplot(reshape2::melt(df_psa_input, variable.name = "Parameter"), 
       aes(x = value)) +
  facet_wrap(~Parameter, scales = "free") +
  geom_histogram(aes(y = ..density..)) +
  ylab("") +
  theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()) 

### 05.05 Run PSA  -------------------------------------------------------------
# Initialize data.frames with PSA output 
# data.frame of costs
df_c <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_c) <- v_names_str

# data.frame of effectiveness
df_e <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_e) <- v_names_str

# Conduct probabilistic sensitivity analysis
# Run Markov model on each parameter set of PSA input dataset
n_time_init_psa_series <- Sys.time()

for (i in 1:n_sim) { # i <- 1
  psa_row <- df_psa_input[i, , drop = FALSE]
  psa_updates <- data.frame(
    name  = names(psa_row),
    value = as.numeric(psa_row[1, ]),
    stringsAsFactors = FALSE
  )
  l_psa_input <- update_param_list(l_params_all, psa_updates)
  
  # Outcomes
  l_out_ce_temp  <- calculate_ce_out(l_psa_input)
  df_c[i, ]  <- l_out_ce_temp$Cost  
  df_e[i, ]  <- l_out_ce_temp$Effect
  
  # Display simulation progress
  if (i/(n_sim/100) == round(i/(n_sim/100), 0)) { # display progress every 5%
    cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}

n_time_end_psa_series <- Sys.time()
n_time_total_psa_series <- n_time_end_psa_series - n_time_init_psa_series

print(paste0("PSA with ", scales::comma(n_sim), 
             " simulations run in series in ", 
             round(n_time_total_psa_series, 2), " ", 
             units(n_time_total_psa_series)))

# ******************************************************************************
# 06 Visualize PSA results for CEA ---------------------------------------------
# ******************************************************************************

### 06.01 Create PSA object  ---------------------------------------------------
l_psa <- dampack::make_psa_obj(
  cost          = df_c, 
  effectiveness = df_e, 
  parameters    = df_psa_input, 
  strategies    = v_names_str
)

l_psa$strategies <- v_names_str
colnames(l_psa$effectiveness) <- v_names_str
colnames(l_psa$cost) <- v_names_str

# Vector with willingness-to-pay (WTP) thresholds
v_wtp <- seq(0, 200000, by = 5000)

### 06.02 Cost-Effectiveness Scatter plot  -------------------------------------
txtsize <- 13

gg_scattter <- plot_psa(l_psa, txtsize = txtsize) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Cost (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")

gg_scattter

### 06.03 Incremental cost-effectiveness ratios (ICERs)  -----------------------
# Compute expected costs and effects for each strategy from the PSA
df_out_ce_psa <- summary(l_psa)

df_cea_psa <- dampack::calculate_icers(
  cost       = df_out_ce_psa$meanCost, 
  effect     = df_out_ce_psa$meanEffect,
  strategies = df_out_ce_psa$Strategy
)

df_cea_psa

### 06.04 Format CEA table  ----------------------------------------------------
# CEA table in proper format 
table_cea <- format_table_cea(df_cea_psa)
table_cea

### 06.05 Plot cost-effectiveness frontier  ------------------------------------
plot_icers(df_cea_psa, label = "all", txtsize = txtsize) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.3))

### 06.06 Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) --
ceac_obj <- dampack::ceac(wtp = v_wtp, psa = l_psa)

# Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)

# CEAC & CEAF plot
gg_ceac <- plot_ceac(ceac_obj, txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.8, 0.48))

gg_ceac

### 06.07 Expected Loss Curves (ELCs)  -----------------------------------------
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj

# ELC plot
gg_elc <- plot_exp_loss(elc_obj, log_y = FALSE, 
                        txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14,
                        col = "full") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Expected Loss (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  theme(legend.position = c(0.4, 0.7))

gg_elc

### 06.08 Expected value of perfect information (EVPI)  ------------------------
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)

# EVPI plot
gg_evpi <- plot_evpi(evpi, effect_units = "QALY", 
                     txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  scale_y_continuous("EVPI (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000)

gg_evpi
