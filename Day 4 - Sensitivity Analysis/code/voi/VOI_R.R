# *****************************************************************************
#
#
# Purpose: Value of Information Analysis using Regression Metamodeling
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
# This code performs Value of Information (VOI) analysis using regression
# metamodeling techniques including:
# - Expected Value of Perfect Information (EVPI)
# - Expected Value of Partial Perfect Information (EVPPI)
# - Expected Value of Sample Information (EVSI)
# - Expected Net Benefit of Sampling (ENBS)
# - Optimal Sample Size (OSS) calculation
#
# *****************************************************************************

# ******************************************************************************
# 01 Setup ---------------------------------------------------------------------
# ******************************************************************************

### 01.01 Clear environment  ---------------------------------------------------
rm(list = ls())

### 01.02 Load packages  -------------------------------------------------------
# Install pacman if not present
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# Load pacman
library(pacman)

# Load (install if needed) CRAN packages
p_load(
  matrixStats,  # Matrix summary statistics
  ggplot2,      # Advanced plotting
  scales,       # Scale functions (dollar labels)
  grid,         # Grid graphics
  dplyr,        # Data manipulation
  reshape2,     # Data reshaping
  mgcv          # Generalized additive models (for fitting splines)
)

### 01.03 Load VOI functions  --------------------------------------------------
source("code/voi/VOI_Functions.R")
source("code/voi/GA_functions.R")

# ******************************************************************************
# 02 Load and prepare data -----------------------------------------------------
# ******************************************************************************

### 02.01 Load PSA simulation file  --------------------------------------------
# Read the `.csv` simulation file into `R`
df_psa <- read.csv("data/PSA.csv", header = TRUE)[, -1]
n_sim  <- nrow(df_psa)

# Display first five observations of the data frame
head(df_psa)

### 02.02 Create Net Monetary Benefit (NMB) matrix  ----------------------------
# Extract NMB columns (WTP = $50,000/QALY)
df_nmb <- df_psa[, 5:7]
head(df_nmb)

# Number of strategies
n_strategies <- ncol(df_nmb)
n_strategies

# Assign name of strategies
strategies <- c("Strategy A", "Strategy B", "Strategy C")
colnames(df_nmb) <- strategies
head(df_nmb)

### 02.03 Visualize NMB distributions  -----------------------------------------
# Format data frame suitably for plotting
df_nmb_long <- reshape2::melt(df_nmb, 
                              variable.name = "Strategy", 
                              value.name = "NMB")

# Plot NMB for different strategies
txtsize <- 16

ggplot(df_nmb_long, aes(x = NMB/1000)) +
  geom_histogram(aes(y = ..density..), col = "black", fill = "gray") +
  geom_density(color = "red") +
  facet_wrap(~ Strategy, scales = "free_y") +
  xlab("Net Monetary Benefit (NMB) per thousand $") +
  scale_x_continuous(n.breaks = 5) + 
  scale_y_continuous("", breaks = NULL, labels = NULL) + 
  theme_bw(base_size = txtsize)

# ******************************************************************************
# 03 Incremental Net Monetary Benefit (INMB) -----------------------------------
# ******************************************************************************

### 03.01 Calculate INMB  ------------------------------------------------------
# Calculate INMB of B vs A
df_inmb <- data.frame(
  Simulation = 1:n_sim,
  `Strategy B vs Strategy A` = df_nmb$`Strategy B` - df_nmb$`Strategy A`
) 

### 03.02 Visualize INMB distribution  -----------------------------------------
# Format data frame suitably for plotting
df_inmb_long <- reshape2::melt(df_inmb, 
                               id.vars = "Simulation", 
                               variable.name = "Comparison", 
                               value.name = "INMB")

# Plot INMB
ggplot(df_inmb_long, aes(x = INMB/1000)) +
  geom_histogram(aes(y = ..density..), col = "black", fill = "gray") +
  geom_density(color = "red") +
  geom_vline(xintercept = 0, col = 4, size = 1.5, linetype = "dashed") +
  facet_wrap(~ Comparison, scales = "free_y") +
  xlab("Incremental Net Monetary Benefit (INMB) in thousand $") +
  scale_x_continuous(n.breaks = 8, limits = c(-100, 100)) + 
  scale_y_continuous("", breaks = NULL, labels = NULL) + 
  theme_bw(base_size = txtsize)

# ******************************************************************************
# 04 Expected Value of Perfect Information (EVPI) ------------------------------
# ******************************************************************************

### 04.01 Find optimal strategy  -----------------------------------------------
# Find optimal strategy (d*) based on the highest expected NMB
d_star <- which.max(colMeans(df_nmb))
d_star

### 04.02 Compute Loss matrix  -------------------------------------------------
# Method 1: Using loop
m_loss <- matrix(0, n_sim, n_strategies)
for (d in 1:n_strategies) { # d <- 1
  m_loss[, d] <- df_nmb[, d] - df_nmb[, d_star]
}
head(m_loss)

# Method 2: Without iterating (much faster!)
m_loss <- as.matrix(df_nmb - df_nmb[, d_star])
head(m_loss)

### 04.03 Calculate EVPI  ------------------------------------------------------
# Find maximum loss overall strategies at each state of the world (i.e., PSA sample)
v_max_loss_i <- rowMaxs(m_loss)
head(v_max_loss_i)

# Average across all states of the world
evpi <- mean(v_max_loss_i)
evpi

# ******************************************************************************
# 05 Expected Value of Partial Perfect Information (EVPPI) ---------------------
# ******************************************************************************

### 05.01 Define parameters  ---------------------------------------------------
v_names_params <- c("Mean No. Visits (A)", 
                    "Mean No. Visits (B)",
                    "Prob. Failing (A)", 
                    "Prob. Failing (B)")

# Matrix with parameters
df_params <- df_psa[, 1:4]
colnames(df_params) <- v_names_params
head(df_params)

# Number and names of parameters
n_params <- ncol(df_params)
n_params

### 05.02 Visualize parameter distributions  -----------------------------------
# Format data suitably for plotting
df_params_long <- reshape2::melt(df_params, variable.name = "Parameter")
head(df_params_long)

# Make parameter names as factors (helps with plotting formatting)
df_params_long$Parameter <- factor(df_params_long$Parameter, 
                                   levels = v_names_params, 
                                   labels = v_names_params)

# Facet plot of parameter distributions
ggplot(df_params_long, aes(x = value)) + 
  geom_histogram(aes(y = ..density..), col = "black", fill = "gray") +
  geom_density(color = "red") +
  facet_wrap(~ Parameter, scales = "free") +
  scale_x_continuous("", n.breaks = 5) + 
  scale_y_continuous("", breaks = NULL, labels = NULL) + 
  theme_bw(base_size = txtsize)

### 05.03 Construct Spline metamodel  ------------------------------------------
# Initialize EVPPI vector 
v_evppi_splines <- matrix(0, n_params)
lmm1 <- vector("list", n_params)
lmm2 <- vector("list", n_params)
lmm3 <- vector("list", n_params)

for (p in 1:n_params) { # p <- 1
  print(paste("Computing EVPPI of parameter", v_names_params[p]))
  
  # Estimate Splines
  lmm1[[p]] <- gam(m_loss[, 1] ~ s(df_params[, p]))
  lmm2[[p]] <- gam(m_loss[, 2] ~ s(df_params[, p]))
  lmm3[[p]] <- gam(m_loss[, 3] ~ s(df_params[, p]))
  
  # Predict Loss using Splines
  m_Lhat_splines <- cbind(lmm1[[p]]$fitted, lmm2[[p]]$fitted, lmm3[[p]]$fitted)
  
  # Compute EVPPI
  v_evppi_splines[p] <- mean(rowMaxs(m_Lhat_splines))
}

### 05.04 Visualize EVPPI results  ---------------------------------------------
# Create data frame for plotting
df_evppi_splines <- data.frame(Parameter = v_names_params, 
                               EVPPI = v_evppi_splines)

df_evppi_splines$Parameter <- factor((df_evppi_splines$Parameter), 
                                     levels = v_names_params[order(df_evppi_splines$EVPPI, 
                                                                   decreasing = TRUE)])

# Plot EVPPI
ggplot(data = df_evppi_splines, aes(x = Parameter, y = EVPPI)) +
  geom_bar(stat = "identity") +
  ylab("EVPPI ($)") +
  scale_y_continuous(n.breaks = 8, labels = comma) +
  theme_bw(base_size = 14)

# ******************************************************************************
# 06 Expected Value of Sample Information (EVSI) -------------------------------
# ******************************************************************************

### 06.01 Define effective (prior) sample size  --------------------------------
n0 <- c(10, # MeanNumVisitsA
        10, # MeanNumVisitsB
        10, # ProbFailA
        10) # ProbFailB

### 06.02 Define future study sample sizes  ------------------------------------
n <- c(0, 1, 5, 10, seq(20, 200, by = 20))
n_samples <- length(n)

### 06.03 Calculate EVSI for each parameter individually  ----------------------
# Initialize EVSI matrix for each parameter
df_evsi <- data.frame(N = n, 
                      matrix(0, nrow = n_samples, ncol = n_params))

# Name columns of EVSI matrix with parameter names
colnames(df_evsi)[-1] <- v_names_params

# Compute EVSI for all parameters separately
for (p in 1:n_params) { # p <- 1
  print(paste("Computing EVSI of parameter", v_names_params[p]))
  
  # Update loss based on Gaussian approximation for each sample of interest
  for (nSamp in 1:n_samples) { # nSamp <- 10
    Ltilde1 <- predict.ga(lmm1[[p]], n = n[nSamp], n0 = n0[p])
    Ltilde2 <- predict.ga(lmm2[[p]], n = n[nSamp], n0 = n0[p])
    Ltilde3 <- predict.ga(lmm3[[p]], n = n[nSamp], n0 = n0[p])
    
    # Combine losses into one matrix
    m_Ltilde <- cbind(Ltilde1, Ltilde2, Ltilde3)
    
    # Apply EVSI equation
    df_evsi[nSamp, p + 1] <- mean(rowMaxs(m_Ltilde))
  }
}

### 06.04 Visualize EVSI results  ----------------------------------------------
# Create EVSI data frame for plotting in decreasing order of EVPPI
df_evsi_long <- reshape2::melt(df_evsi[1:21, ], 
                               id.vars = "N", 
                               variable.name = "Parameter", 
                               value.name = "evsi")

df_evsi_long$Parameter <- factor((df_evsi_long$Parameter), 
                                 levels = v_names_params[order(df_evppi_splines$EVPPI, 
                                                               decreasing = TRUE)])

# Plot EVSI
ggplot(df_evsi_long, aes(x = N, y = evsi)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Parameter) +
  ggtitle("Expected Value of Sample Information (EVSI)") +
  scale_x_continuous("Sample size (n)", n.breaks = 5) + 
  scale_y_continuous("", n.breaks = 6, labels = dollar) + 
  theme_bw(base_size = txtsize)

# Plot EVSI with EVPPI reference line
ggplot(df_evsi_long, aes(x = N, y = evsi)) +
  geom_line(aes(linetype = "EVSI")) +
  geom_point() +
  facet_wrap(~ Parameter) +
  geom_hline(aes(yintercept = EVPPI, 
                 linetype = "EVPPI"), 
             data = df_evppi_splines) +
  scale_linetype_manual(name = "", 
                        values = c("EVSI" = "solid", "EVPPI" = "dashed")) +
  scale_x_continuous("Sample size (n)", n.breaks = 5) + 
  scale_y_continuous("", n.breaks = 6, labels = dollar) + 
  theme_bw(base_size = txtsize) +
  theme(legend.position = c(0.8, 0.7))

# ******************************************************************************
# 07 EVSI for combination of parameters ----------------------------------------
# ******************************************************************************

### 07.01 Observational study design  ------------------------------------------
#### Define parameters and sample sizes
v_sel_params_obs <- c(1, 2)

# Vector with samples to evaluate EVSI for an Observational design
v_n_obs <- c(0, 1, 5, 10, seq(20, 200, by = 20), 300, 400, 500, 600, 700, 800)
n_obs_samples <- length(v_n_obs)

# Initialize EVSI matrix for a combination of parameters
df_evsi_obs <- data.frame(Study = "Observational", 
                          N = v_n_obs, 
                          EVSI = matrix(0, nrow = n_obs_samples, ncol = 1))

#### Estimate linear metamodel of two parameters
lmm1_obs <- gam(m_loss[, 1] ~ s(df_params[, v_sel_params_obs[1]]) + 
                  s(df_params[, v_sel_params_obs[2]]) + 
                  ti(df_params[, v_sel_params_obs[1]], 
                     df_params[, v_sel_params_obs[2]]))

lmm2_obs <- gam(m_loss[, 2] ~ s(df_params[, v_sel_params_obs[1]]) + 
                  s(df_params[, v_sel_params_obs[2]]) + 
                  ti(df_params[, v_sel_params_obs[1]], 
                     df_params[, v_sel_params_obs[2]]))

lmm3_obs <- gam(m_loss[, 3] ~ s(df_params[, v_sel_params_obs[1]]) + 
                  s(df_params[, v_sel_params_obs[2]]) + 
                  ti(df_params[, v_sel_params_obs[1]], 
                     df_params[, v_sel_params_obs[2]]))

# Predict Loss using Splines
m_Lhat_obs_splines <- cbind(lmm1_obs$fitted, lmm2_obs$fitted, lmm3_obs$fitted)

# Compute EVPPI
evppi_obs <- mean(rowMaxs(m_Lhat_obs_splines))
evppi_obs          

#### Compute EVSI for groups of parameters
for (nSamp in 1:n_obs_samples) {
  Ltilde1_obs <- predict.ga(lmm1_obs, n = v_n_obs[nSamp], 
                            n0 = n0[v_sel_params_obs])
  Ltilde2_obs <- predict.ga(lmm2_obs, n = v_n_obs[nSamp], 
                            n0 = n0[v_sel_params_obs])
  Ltilde3_obs <- predict.ga(lmm3_obs, n = v_n_obs[nSamp], 
                            n0 = n0[v_sel_params_obs])
  
  # Combine losses into one matrix
  m_Ltilde_obs <- cbind(Ltilde1_obs, Ltilde2_obs, Ltilde3_obs)
  
  # Apply EVSI equation
  df_evsi_obs$EVSI[nSamp] <- mean(rowMaxs(m_Ltilde_obs))
}

### 07.02 Randomized Controlled Trial (RCT) design  ----------------------------
#### Define parameters and sample sizes
v_sel_params_rct <- c(3, 4)

# Vector with samples to evaluate EVSI for a RCT
v_n_rct <- c(0, 1, 5, 10, seq(20, 200, by = 20))
n_rct_samples <- length(v_n_rct)

# Initialize EVSI matrix for a combination of parameters
df_evsi_rct <- data.frame(Study = "RCT",
                          N = v_n_rct, 
                          EVSI = matrix(0, nrow = n_rct_samples, ncol = 1))

#### Estimate linear metamodel of two parameters
lmm1_rct <- gam(m_loss[, 1] ~ s(df_params[, v_sel_params_rct[1]]) + 
                  s(df_params[, v_sel_params_rct[2]]) + 
                  ti(df_params[, v_sel_params_rct[1]], 
                     df_params[, v_sel_params_rct[2]]))

lmm2_rct <- gam(m_loss[, 2] ~ s(df_params[, v_sel_params_rct[1]]) + 
                  s(df_params[, v_sel_params_rct[2]]) + 
                  ti(df_params[, v_sel_params_rct[1]], 
                     df_params[, v_sel_params_rct[2]]))

lmm3_rct <- gam(m_loss[, 3] ~ s(df_params[, v_sel_params_rct[1]]) + 
                  s(df_params[, v_sel_params_rct[2]]) + 
                  ti(df_params[, v_sel_params_rct[1]], 
                     df_params[, v_sel_params_rct[2]]))

# Predict Loss using Splines
m_Lhat_rct_splines <- cbind(lmm1_rct$fitted, lmm2_rct$fitted, lmm3_rct$fitted)

# Compute EVPPI
evppi_rct <- mean(rowMaxs(m_Lhat_rct_splines))
evppi_rct          

#### Compute EVSI over different sample sizes
for (nSamp in 1:n_rct_samples) {
  Ltilde1_rct <- predict.ga(lmm1_rct, n = v_n_rct[nSamp], n0 = n0[v_sel_params_rct])
  Ltilde2_rct <- predict.ga(lmm2_rct, n = v_n_rct[nSamp], n0 = n0[v_sel_params_rct])
  Ltilde3_rct <- predict.ga(lmm3_rct, n = v_n_rct[nSamp], n0 = n0[v_sel_params_rct])
  
  # Combine losses into one matrix
  m_Ltilde_rct <- cbind(Ltilde1_rct, Ltilde2_rct, Ltilde3_rct)
  
  # Apply EVSI equation
  df_evsi_rct$EVSI[nSamp] <- mean(rowMaxs(m_Ltilde_rct))
}

### 07.03 Visualize EVSI for both study designs  -------------------------------
# Combine both study designs
df_evppi_combo <- data.frame(Study = c("Observational", "RCT"), 
                             EVPPI = c(evppi_obs, evppi_rct))

df_evsi_combo <- bind_rows(df_evsi_obs, df_evsi_rct)

# Plot EVSI by study design
ggplot(df_evsi_combo, aes(x = N, y = EVSI)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Study, scales = "free_x") +
  geom_hline(aes(yintercept = EVPPI, 
                 linetype = "EVPPI"), 
             data = df_evppi_combo) +
  scale_linetype_manual(name = "", 
                        values = c("EVSI" = "solid", "EVPPI" = "dashed")) +
  ggtitle("EVSI for different study designs") +
  scale_x_continuous("Sample size (n)", n.breaks = 5) + 
  scale_y_continuous("", n.breaks = 8, labels = dollar) + 
  theme_bw(base_size = txtsize) +
  theme(legend.position = c(0.2, 0.9))

# ******************************************************************************
# 08 Expected Net Benefit of Sampling (ENBS) -----------------------------------
# ******************************************************************************

### 08.01 Define population values  --------------------------------------------
# Discount rate
disc <- c(0.03)

# Technology lifetime
LT   <- 10
v_time <- seq(0, LT)

# Annual Number of Individuals to Be Treated
# Present prevalence (in millions)
prev  <- 0.010

# Annual Incidence (in millions)
incid <- 147*1e-6

# Total population affected by technology calculated with `TotPop` function
tot_pop <- TotPop(v_time, prev, incid, disc) 

### 08.02 Calculate Population EVSI  -------------------------------------------
# Observational study
df_pop_evsi_obs <- df_evsi_obs
df_pop_evsi_obs$popEVSI <- df_pop_evsi_obs$EVSI * tot_pop

# RCT
df_pop_evsi_rct <- df_evsi_rct
df_pop_evsi_rct$popEVSI <- df_pop_evsi_rct$EVSI * tot_pop

### 08.03 Calculate cost of research  ------------------------------------------
# Observational study
v_cost_res_obs <- CostRes(
  fixed.cost = 10000e-6,
  samp.size = v_n_obs,
  cost.per.patient = 500e-6,  # In Million $
  INMB = 0,
  clin.trial = FALSE
)

# Data frame with cost of trial in Millions
df_cost_obs <- data.frame(N = v_n_obs, CS = v_cost_res_obs)

# RCT
v_cost_res_rct <- CostRes(
  fixed.cost = 8000000e-6,
  samp.size = v_n_rct,
  cost.per.patient = 8500e-6,  # In Million $
  INMB = 0,
  clin.trial = TRUE
) 

# Data frame with cost of trial in Millions
df_cost_rct <- data.frame(N = v_n_rct, CS = v_cost_res_rct)

### 08.04 Calculate ENBS  ------------------------------------------------------
# Create ENBS data frame
df_enbs_obs <- merge(df_pop_evsi_obs, df_cost_obs, by = "N")
df_enbs_rct <- merge(df_pop_evsi_rct, df_cost_rct, by = "N")

# Compute ENBS 
df_enbs_obs$ENBS <- df_enbs_obs$popEVSI - df_enbs_obs$CS
df_enbs_rct$ENBS <- df_enbs_rct$popEVSI - df_enbs_rct$CS

# Compute OSS (n*)
df_enbs_obs$nstar <- df_enbs_obs$N[which.max(df_enbs_obs$ENBS)]
df_enbs_rct$nstar <- df_enbs_rct$N[which.max(df_enbs_rct$ENBS)]

# Append data frames
df_enbs_all <- bind_rows(df_enbs_obs, df_enbs_rct)

### 08.05 Calculate Optimal Sample Size (OSS)  ---------------------------------
df_oss <- summarise(
  group_by(df_enbs_all, Study),
  MaxENBS = max(ENBS),
  Nstar   = N[which.max(ENBS)]
)
df_oss

### 08.06 Visualize ENBS results  ----------------------------------------------
# Create suitable data frames for plotting
df_enbs_obs_long <- reshape2::melt(df_enbs_obs[, -3], 
                                   id.vars = c("Study", "N", "nstar"), 
                                   value.name = "Million")

df_enbs_rct_long <- reshape2::melt(df_enbs_rct[, -3], 
                                   id.vars = c("Study", "N", "nstar"), 
                                   value.name = "Million")

# Append data frames for plotting
df_enbs_ll_long <- bind_rows(df_enbs_obs_long, df_enbs_rct_long)

levels(df_enbs_ll_long$Study) <- c(
  paste("Observational; n* = ", comma(df_oss$Nstar[1]), sep = ""), 
  paste("RCT; n* = ", comma(df_oss$Nstar[2]), sep = "")
)

# Plot ENBS, popEVSI, and cost of research
ggplot(df_enbs_ll_long, aes(x = N, y = Million, 
                            colour = variable, group = variable)) + 
  facet_wrap(~ Study, scales = "free_x") +
  geom_hline(aes(yintercept = 0), size = 0.7, 
             linetype = 2, colour = "gray") + 
  geom_vline(aes(xintercept = nstar), size = 0.7, 
             linetype = 2, colour = "gray") + 
  geom_point() +
  geom_line() +
  scale_x_continuous("Sample size (N)", n.breaks = 6, labels = comma) +
  scale_y_continuous("Value (Million $)", n.breaks = 6, labels = comma, 
                     limits = c(0, 40)) +
  scale_colour_hue("Study design ", l = 50,
                   labels = c("popEVSI(n) ", 
                              "Cost of Research(n) ", 
                              "ENBS(n) ")) +
  theme_bw(base_size = txtsize) +
  theme(legend.position = "bottom",
        panel.margin = unit(2, "lines"))

# ******************************************************************************
# 09 Summary: Optimal Sample Size (OSS) ----------------------------------------
# ******************************************************************************

### 09.01 Display OSS results  -------------------------------------------------
df_oss
