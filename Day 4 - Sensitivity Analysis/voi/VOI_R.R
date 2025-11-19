# title: "Value of Information Analysis using Regression Metamodeling"

## Clean everything from the workspace
rm(list = ls())

## Install Packages (If you haven't done so yet!)
# install.packages("gdata")
# install.packages("xlsx")
# install.packages("ggplot2")
# install.packages("reshape2")
# install.packages("scales")
# install.packages("grid")
# install.packages("gridExtra")
# install.packages("dplyr")
# install.packages("matrixStats")
library(matrixStats)
library(ggplot2)
library(scales)  # For dollar labels
library(grid)
library(dplyr)
library(reshape2)
library(mgcv) # For fitting splines

#### Load VOI Functions ####
source("R-labs/value-of-information/VOI_Functions.R")
source("R-labs/value-of-information/GA_functions.R")

#### Simple Example: WTP = $50,000/QALY ####
## Load simulation file
# Read the `.csv` simulation file into `R`.
df_psa <- read.csv("R-labs/value-of-information/data/PSA.csv", header = TRUE)[, -1]
n_sim  <- nrow(df_psa)

#Display first five observations of the data frame using the command `head`
head(df_psa)

### Net Monetary Benefit (NMB) ####
# Create NMB matrix
df_nmb <- df_psa[, 5:7]
head(df_nmb)

# Number of Strategies
n_strategies <- ncol(df_nmb)
n_strategies

# Assign name of strategies
strategies <- c("Strategy A", "Strategy B", "Strategy C")
colnames(df_nmb) <- strategies
head(df_nmb)

## Format data frame suitably for plotting
df_nmb_long <- reshape2::melt(df_nmb, 
                              variable.name = "Strategy", 
                              value.name = "NMB")

## Plot NMB for different strategies
txtsize <- 16
# Faceted plot by Strategy
ggplot(df_nmb_long, aes(x = NMB/1000)) +
  geom_histogram(aes(y =..density..), col="black", fill = "gray") +
  geom_density(color = "red") +
  facet_wrap(~ Strategy, scales = "free_y") +
  xlab("Net Monetary Benefit (NMB) per thousand $") +
  scale_x_continuous(n.breaks = 5) + 
  scale_y_continuous("", breaks = NULL, labels = NULL) + 
  theme_bw(base_size = txtsize)

#### Incremental NMB (INMB) ####
# Calculate INMB of B vs A
# Only B vs A but we could have plotted all combinations
df_inmb <- data.frame(Simulation = 1:n_sim,
                      `Strategy B vs Strategy A` = df_nmb$`Strategy B` - df_nmb$`Strategy A`) 

## Format data frame suitably for plotting
df_inmb_long <- reshape2::melt(df_inmb, id.vars = "Simulation", 
                               variable.name = "Comparison", 
                               value.name = "INMB")
## Plot INMB
ggplot(df_inmb_long, aes(x = INMB/1000)) +
  geom_histogram(aes(y = ..density..), col = "black", fill = "gray") +
  geom_density(color = "red") +
  geom_vline(xintercept = 0, col = 4, size = 1.5, linetype = "dashed") +
  facet_wrap(~ Comparison, scales = "free_y") +
  xlab("Incremental Net Monetary Benefit (INMB) in thousand $") +
  scale_x_continuous(n.breaks = 8, limits = c(-100, 100)) + 
  scale_y_continuous("", breaks = NULL, labels = NULL) + 
  theme_bw(base_size = txtsize)

#### Loss Matrix ####  
# Find optimal strategy (d*) based on the highest expected NMB
d_star <- which.max(colMeans(df_nmb))
d_star

## Compute Loss matrix iterating over all strategies
# Initialize loss matrix of dimension: number of simulation by number of strategies
m_loss <- matrix(0, n_sim, n_strategies)
for (d in 1:n_strategies) { # d <- 1
  m_loss[, d] <- df_nmb[, d] - df_nmb[, d_star]
}
head(m_loss)

# Or without iterating (much faster!)
m_loss <- as.matrix(df_nmb - df_nmb[, d_star])
head(m_loss)

# EVPI ----
## Find maximum loss overall strategies at each state of the world 
## (i.e., PSA sample)
v_max_loss_i <- rowMaxs(m_loss)
head(v_max_loss_i)
## Average across all states of the world
evpi <- mean(v_max_loss_i)
evpi

# EVPPI ----
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

### Histogram of parameters
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

## Construct Spline metamodel ----
### Splines
## Initialize EVPPI vector 
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

## Plotting EVPPI using order of polynomial
df_evppi_splines <- data.frame(Parameter = v_names_params, 
                               EVPPI = v_evppi_splines)
df_evppi_splines$Parameter <- factor((df_evppi_splines$Parameter), 
                              levels = v_names_params[order(df_evppi_splines$EVPPI, 
                                                            decreasing = TRUE)])

# Plot EVPPI using ggplot2 package
ggplot(data = df_evppi_splines, aes(x = Parameter, y = EVPPI)) +
  geom_bar(stat = "identity") +
  ylab("EVPPI ($)") +
  scale_y_continuous(n.breaks = 8, labels = comma) +
  theme_bw(base_size = 14)

# EVSI ----
## Effective (prior) Sample size ----
n0 <- c(10, # MeanNumVisitsA
        10, # MeanNumVisitsB
        10, # ProbFailA
        10) # ProbFailB

## Future study sample sizes ----
n <- c(0, 1, 5, 10, seq(20, 200, by = 20))
n_samples <- length(n)

## Each parameter individually (only assuming linear relationship) ----
# Initialize EVSI matrix for each parameters
df_evsi <- data.frame(N = n, 
                   matrix(0, nrow = n_samples, ncol = n_params))

# Name columns of EVPSI matrix with parameter names
colnames(df_evsi)[-1] <- v_names_params

### Compute EVSI for all parameters separately ----
for (p in 1:n_params) { # p <- 1
  print(paste("Computing EVSI of parameter", v_names_params[p]))
    # Update loss based on Gaussian approximation for each sample of interest
    for (nSamp in 1:n_samples) { # nSamp <- 10
      Ltilde1 <- predict.ga(lmm1[[p]], n = n[nSamp], n0 = n0[p])
      Ltilde2 <- predict.ga(lmm2[[p]], n = n[nSamp], n0 = n0[p])
      Ltilde3 <- predict.ga(lmm3[[p]], n = n[nSamp], n0 = n0[p])
      ## Combine losses into one matrix
      m_Ltilde <- cbind(Ltilde1, Ltilde2, Ltilde3)
      ### Apply EVSI equation
      df_evsi[nSamp, p + 1] <- mean(rowMaxs(m_Ltilde))
    }
}

### Plotting EVSI ----
# Create EVSI data frame for plotting in decreasing order of EVPPI
df_evsi_long <- reshape2::melt(df_evsi[1:21,], 
                               id.vars = "N", 
                               variable.name = "Parameter", 
                               value.name = "evsi")
df_evsi_long$Parameter <- factor((df_evsi_long$Parameter), 
                                 levels = v_names_params[order(df_evppi_splines$EVPPI, 
                                                               decreasing = TRUE)])

# Plot evsi using ggplot2 package
ggplot(df_evsi_long, aes(x = N, y = evsi)) +  # colour = Parameter
  geom_line() +
  geom_point() +
  facet_wrap(~ Parameter) +  # scales = "free_y"
  ggtitle("Expected Value of Sample Information (EVSI)") +
  scale_x_continuous("Sample size (n)", n.breaks = 5) + 
  scale_y_continuous("", n.breaks = 6, labels = dollar) + 
  theme_bw(base_size = txtsize)
# ggsave("Figs/Toy_evsi.pdf", width = 8, height = 6)
# ggsave("Figs/Toy_evsi.png", width = 8, height = 6)

# Adding EVPPI 
ggplot(df_evsi_long, aes(x = N, y = evsi)) +  # colour = Parameter
  geom_line(aes(linetype = "EVSI")) +
  geom_point() +
  facet_wrap(~ Parameter) +  # scales = "free_y"
  geom_hline(aes(yintercept = EVPPI, 
                 linetype = "EVPPI"), 
             data = df_evppi_splines) +
  scale_linetype_manual(name = "", 
                        values = c("EVSI" = "solid", "EVPPI" = "dashed")) +
  scale_x_continuous("Sample size (n)", n.breaks = 5) + 
  scale_y_continuous("", n.breaks = 6, labels = dollar) + 
  theme_bw(base_size = txtsize) +
  theme(legend.position = c(0.8, 0.7))
# ggsave("Figs/Toy_EVPSI_EVPPI.pdf", width = 8, height = 6)
# ggsave("Figs/Toy_EVPSI_EVPPI.png", width = 8, height = 6)

## Combination of parameters ----
### Assuming an observational study ----
v_sel_params_obs <- c(1, 2)
# Vector with samples to evaluate EVPSI for an Observational design
v_n_obs <- c(0, 1, 5, 10, seq(20, 200, by = 20), 300, 400, 500, 600, 700, 800) #seq(0, 1000, by = 20)
n_obs_samples <- length(v_n_obs)
# Initialize EVPSI matrix for a combination of parameters
df_evsi_obs <- data.frame(Study = "Observational", 
                          N = v_n_obs, 
                          EVSI = matrix(0, nrow = n_obs_samples, ncol = 1))

#### Estimate linear metamodel of two parameters ----
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

#### Compute EVSI for groups of parameters ----
for (nSamp in 1:n_obs_samples) {
  Ltilde1_obs <- predict.ga(lmm1_obs, n = v_n_obs[nSamp], 
                            n0 = n0[v_sel_params_obs])
  Ltilde2_obs <- predict.ga(lmm2_obs, n = v_n_obs[nSamp], 
                            n0 = n0[v_sel_params_obs])
  Ltilde3_obs <- predict.ga(lmm3_obs, n = v_n_obs[nSamp], 
                            n0 = n0[v_sel_params_obs])
  ## Combine losses into one matrix
  m_Ltilde_obs <- cbind(Ltilde1_obs, Ltilde2_obs, Ltilde3_obs)
  ### Apply EVSI equation
  df_evsi_obs$EVSI[nSamp] <- mean(rowMaxs(m_Ltilde_obs))
}

### Assuming an RCT ----
v_sel_params_rct <- c(3, 4)
# Vector with samples to evaluate EVPSI for a RCT
v_n_rct <- c(0, 1, 5, 10, seq(20, 200, by = 20))
n_rct_samples <- length(v_n_rct)
# Initialize EVPSI matrix for a combination of parameters
df_evsi_rct <- data.frame(Study = "RCT",
                          N = v_n_rct, 
                          EVSI = matrix(0, nrow = n_rct_samples, ncol = 1))

#### Estimate linear metamodel of two parameters ----
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

#### Compute EVSI over different sample sizes ----
for (nSamp in 1:n_rct_samples) {
  Ltilde1_rct <- predict.ga(lmm1_rct, n = v_n_rct[nSamp], n0 = n0[v_sel_params_rct])
  Ltilde2_rct <- predict.ga(lmm2_rct, n = v_n_rct[nSamp], n0 = n0[v_sel_params_rct])
  Ltilde3_rct <- predict.ga(lmm3_rct, n = v_n_rct[nSamp], n0 = n0[v_sel_params_rct])
  ## Combine losses into one matrix
  m_Ltilde_rct <- cbind(Ltilde1_rct, Ltilde2_rct, Ltilde3_rct)
  ### Apply EVSI equation
  df_evsi_rct$EVSI[nSamp] <- mean(rowMaxs(m_Ltilde_rct))
}


## Plot EVSI for both study designs
# Combine both study designs
df_evppi_combo <- data.frame(Study = c("Observational", "RCT"), 
                             EVPPI = c(evppi_obs, evppi_rct))
df_evsi_combo <- bind_rows(df_evsi_obs,
                           df_evsi_rct)

# Plot EVSI by study design
ggplot(df_evsi_combo, aes(x = N, y = EVSI)) +  # colour = Parameter
  geom_line() +
  geom_point() +
  facet_wrap(~ Study, scales = "free_x") +
  geom_hline(aes(yintercept = EVPPI, 
                 linetype = "EVPPI"), 
             data = df_evppi_combo) +
  scale_linetype_manual(name = "", 
                        values = c("EVSI" = "solid", "EVPPI" = "dashed")) +
  ggtitle("EVPSI for different study designs") +
  scale_x_continuous("Sample size (n)", n.breaks = 5) + 
  scale_y_continuous("", n.breaks = 8, labels = dollar) + 
  theme_bw(base_size = txtsize) +
  theme(legend.position = c(0.2, 0.9))
# ggsave("Figs/Toy_EVPSI_Studies.pdf", width = 8, height = 6)
# ggsave("Figs/Toy_EVPSI_Studies.png", width = 8, height = 6)

#### ENBS ####
### Population Values
## Discount rate
disc <- c(0.03)
# Technology lifetime
LT   <- 10
v_time <- seq(0, LT)
## Annual Number of Individuals to Be Treated
# Present prevalence
prev     <- 0.010 # In millions(1e6)
# Annual Incidence. 
incid    <- 147*1e-6 # In millions: 0.005*29.376e-3
## Total population afectd by technology calculated with `TotPop` function in Millions
tot_pop <- TotPop(v_time,    # Function
                  prev, 
                  incid, 
                  disc) 
## Population EVPSI
# Observational study
df_pop_evsi_obs <- df_evsi_obs
df_pop_evsi_obs$popEVSI <- df_pop_evsi_obs$EVSI*tot_pop
# RCT
df_pop_evsi_rct <- df_evsi_rct
df_pop_evsi_rct$popEVSI <- df_pop_evsi_rct$EVSI*tot_pop

### Cost of research
## Observational study
v_cost_res_obs <- CostRes(fixed.cost = 10000e-6,
                          samp.size = v_n_obs,  # vector 
                          cost.per.patient = 500e-6, # In Million $
                          INMB = 0,
                          clin.trial = FALSE)
# Data frame with cost of trial in Millions
df_cost_obs <- data.frame(N = v_n_obs, 
                          CS = v_cost_res_obs)
## RCT
v_cost_res_rct <- CostRes(fixed.cost = 8000000e-6,
                        samp.size = v_n_rct,  # vector 
                        cost.per.patient = 8500e-6, # In Million $
                        INMB = 0,
                        clin.trial = TRUE) 
# Data frame with cost of trial in Millions
df_cost_rct <- data.frame(N = v_n_rct, CS = v_cost_res_rct)

### Create ENBS data frame
df_enbs_obs <- merge(df_pop_evsi_obs, df_cost_obs, by = "N")
df_enbs_rct <- merge(df_pop_evsi_rct, df_cost_rct, by = "N")

## Compute ENBS 
df_enbs_obs$ENBS <- df_enbs_obs$popEVSI - df_enbs_obs$CS
df_enbs_rct$ENBS <- df_enbs_rct$popEVSI - df_enbs_rct$CS

## Compute OSS (n*)
df_enbs_obs$nstar <- df_enbs_obs$N[which.max(df_enbs_obs$ENBS)]
df_enbs_rct$nstar <- df_enbs_rct$N[which.max(df_enbs_rct$ENBS)]

# Append data frames
df_enbs_all <- bind_rows(df_enbs_obs,
                         df_enbs_rct)

df_oss <- summarise(group_by(df_enbs_all, Study),
                    MaxENBS = max(ENBS),
                    Nstar   = N[which.max(ENBS)])
df_oss

## Plot ENBS, EVPSI and n*
# Create suitable data frames for plotting
df_enbs_obs_long <- reshape2::melt(df_enbs_obs[, -3], 
                                   id.vars = c("Study", "N", "nstar"), 
                                   value.name = "Million")
df_enbs_rct_long <- reshape2::melt(df_enbs_rct[, -3], 
                                   id.vars = c("Study", "N", "nstar"), 
                                   value.name = "Million")
# Append data frames for plotting
df_enbs_ll_long <- bind_rows(df_enbs_obs_long,
                             df_enbs_rct_long)
levels(df_enbs_ll_long$Study) <- c(paste("Observational; n* = ", 
                                         comma(df_oss$Nstar[1]), sep = ""), 
                                   paste("RCT; n* = ", 
                                         comma(df_oss$Nstar[2]), sep = ""))
ggplot(df_enbs_ll_long, aes(x = N, y = Million, 
                            colour = variable, group = variable)) + 
  facet_wrap(~ Study, scales = "free_x") +
  #geom_segment(data = oss, aes(x = Nstar, y = 0, xend = Nstar, yend = MaxENBS)) + 
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
                   labels = c("popEVPSI(n) ", 
                              "Cost of Research(n) ", 
                              "ENBS(n) ")) +
  theme_bw(base_size = txtsize) +
  theme(legend.position = "bottom",
        panel.margin = unit(2, "lines"))
# ggsave("Figs/Toy_ENBS.pdf", width = 8, height = 6)
# ggsave("Figs/Toy_ENBS.png", width = 8, height = 6)

#### Optimal Sample Size (OSS), n*
df_oss

