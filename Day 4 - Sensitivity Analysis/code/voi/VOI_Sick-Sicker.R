# *****************************************************************************
#
#
# Purpose: Value of Information Analysis for the Sick-Sicker model using 
#          Regression Metamodeling
#
# Authors: 
# This work is developed by the Decision Analysis in R for Technologies in Health 
# (DARTH) workgroup:
#
# - Hawre J. Jalal, MD, PhD
# - Fernando Alarid-Escudero, PhD <falarid@stanford.edu>
# - Eva A. Enns, MS, PhD 
# - M.G. Myriam Hunink, MD, PhD 
# - Eline Krijkamp, PhD 
# - Petros Pechlivanoglou, PhD
# - Alan Yang, MSc
#
# *****************************************************************************
#
# Notes:
#
# This code performs Value of Information (VOI) analysis for the Sick-Sicker
# model using regression metamodeling techniques including:
# - Expected Value of Perfect Information (EVPI)
# - Expected Value of Partial Perfect Information (EVPPI)
# - Expected Value of Sample Information (EVSI)
#
# Presented at: ESMDM 2018 Short Course, June 10th, Leiden, The Netherlands
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
toy <- read.csv("data/psa_sick_sicker.csv", header = TRUE)[, -1]

n.sim <- nrow(toy)

# Display first five observations of the data frame
head(toy)

### 02.02 Create Net Monetary Benefit (NMB) matrix  ----------------------------
# Define willingness-to-pay threshold
wtp <- 120000

# Calculate NMB for each strategy
toy$NMB_NoTrt <- wtp * toy$QALY_NoTrt - toy$Cost_NoTrt
toy$NMB_Trt   <- wtp * toy$QALY_Trt - toy$Cost_Trt

# Extract NMB columns
nmb <- toy[, c("NMB_NoTrt", "NMB_Trt")]
head(nmb)

# Number of strategies
n.strategies <- ncol(nmb)
n.strategies

# Assign name of strategies
strategies <- c("No Trt", "Trt")
colnames(nmb) <- strategies
head(nmb)

### 02.03 Visualize NMB distributions  -----------------------------------------
# Format data frame suitably for plotting
nmb.gg <- melt(nmb,  
               variable.name = "Strategy", 
               value.name = "NMB")

# Plot NMB for different strategies
ggplot(nmb.gg, aes(x = NMB/1000)) +
  geom_histogram(aes(y = ..density..), col = "black", fill = "gray") +
  geom_density(color = "red") +
  facet_wrap(~ Strategy, scales = "free_y") +
  xlab("Net Monetary Benefit (NMB) x10^3") +
  scale_x_continuous(breaks = number_ticks(5), labels = dollar) + 
  scale_y_continuous(breaks = number_ticks(5)) + 
  theme_bw()

# ******************************************************************************
# 03 Incremental Net Monetary Benefit (INMB) -----------------------------------
# ******************************************************************************

### 03.01 Calculate INMB  ------------------------------------------------------
# Calculate INMB of Trt vs. No Trt
inmb <- data.frame(
  Simulation = 1:n.sim,
  `Trt vs. No Trt` = nmb$Trt - nmb$`No Trt`
) 

### 03.02 Visualize INMB distribution  -----------------------------------------
# Format data frame suitably for plotting
inmb.gg <- melt(inmb, 
                id.vars = "Simulation", 
                variable.name = "Comparison", 
                value.name = "INMB")

txtsize <- 16

# Plot INMB
ggplot(inmb.gg, aes(x = INMB/1000)) +
  geom_histogram(aes(y = ..density..), col = "black", fill = "gray") +
  geom_density(color = "red") +
  geom_vline(xintercept = 0, col = 4, size = 1.5, linetype = "dashed") +
  facet_wrap(~ Comparison, scales = "free_y") +
  xlab("Incremental Net Monetary Benefit (INMB) in thousand $") +
  scale_x_continuous(breaks = number_ticks(5), limits = c(-100, 100)) + 
  scale_y_continuous(breaks = number_ticks(5)) + 
  theme_bw(base_size = txtsize)

# ******************************************************************************
# 04 Expected Value of Perfect Information (EVPI) ------------------------------
# ******************************************************************************

### 04.01 Find optimal strategy  -----------------------------------------------
# Find optimal strategy (d*) based on the highest expected NMB
d.star <- which.max(colMeans(nmb))
d.star

### 04.02 Compute Loss matrix  -------------------------------------------------
loss <- as.matrix(nmb - nmb[, d.star])
head(loss)

### 04.03 Calculate EVPI  ------------------------------------------------------
# Find maximum loss overall strategies at each state of the world (i.e., PSA sample)
max.loss.i <- rowMaxs(loss)
head(max.loss.i)

# Average across all states of the world
evpi <- mean(max.loss.i)
evpi

# ******************************************************************************
# 05 Expected Value of Partial Perfect Information (EVPPI) ---------------------
# ******************************************************************************

### 05.01 Define parameters  ---------------------------------------------------
# Matrix with parameters
x <- toy[, c(1:14)]
head(x)

# Number and names of parameters
n.params <- ncol(x)
n.params

names.params <- colnames(x) 
names.params

### 05.02 Visualize parameter distributions  -----------------------------------
# Format data suitably for plotting
params <- melt(x, variable.name = "Parameter")
head(params)

# Make parameter names as factors (helps with plotting formatting)
params$Parameter <- factor(params$Parameter, 
                           levels = names.params, 
                           labels = names.params)

# Facet plot of parameter distributions
ggplot(params, aes(x = value)) + 
  geom_histogram(aes(y = ..density..), col = "black", fill = "gray") +
  geom_density(color = "red") +
  facet_wrap(~ Parameter, scales = "free") +
  scale_x_continuous(breaks = number_ticks(5)) + 
  scale_y_continuous(breaks = number_ticks(5)) + 
  theme_bw(base_size = 14)

### 05.03 Construct Spline metamodel  ------------------------------------------
# Initialize EVPPI vector 
evppi.splines <- matrix(0, n.params)
lmm1 <- vector("list", n.params)
lmm2 <- vector("list", n.params)

for (p in 1:n.params) { # p <- 1
  print(paste("Computing EVPPI of parameter", names.params[p]))
  
  # Estimate Splines
  lmm1[[p]] <- gam(loss[, 1] ~ s(x[, p]))
  lmm2[[p]] <- gam(loss[, 2] ~ s(x[, p]))
  
  # Predict Loss using Splines
  Lhat.splines <- cbind(lmm1[[p]]$fitted, lmm2[[p]]$fitted)
  
  # Compute EVPPI
  evppi.splines[p] <- mean(rowMaxs(Lhat.splines))
}

### 05.04 Visualize EVPPI results  ---------------------------------------------
# Create data frame for plotting
evppi.splines.gg <- data.frame(Parameter = names.params, 
                               EVPPI = evppi.splines)

evppi.splines.gg$Parameter <- factor((evppi.splines.gg$Parameter), 
                                     levels = names.params[order(evppi.splines.gg$EVPPI, 
                                                                 decreasing = TRUE)])

# Plot EVPPI
ggplot(data = evppi.splines.gg, aes(x = Parameter, y = EVPPI)) +
  geom_bar(stat = "identity") +
  ylab("EVPPI ($)") +
  scale_y_continuous(breaks = number_ticks(6), labels = comma) +
  theme_bw(base_size = 14)

# ******************************************************************************
# 06 Expected Value of Sample Information (EVSI) -------------------------------
# ******************************************************************************

### 06.01 Select parameters with positive EVPPI  -------------------------------
sel.params <- c(3, 4, 10, 12, 14)
n.params <- length(sel.params)

### 06.02 Define effective (prior) sample sizes  -------------------------------
# Effective (prior) Sample size
n0 <- numeric(length(sel.params))
n0[1] <- 84 + 800    # p.S1S2 ~ Beta(84, 800)
n0[2] <- 10 + 2000   # p.HD ~ Beta(10, 2000)
n0[3] <- 73.5        # cTrt ~ Gamma(73.5, 163.3) -> likelihood ~ Exponential
n0[4] <- 50          # u.S1 ~ N(0.75, 0.02 / sqrt(50))
n0[5] <- 20          # u.Trt ~ N(0.95, 0.02)

### 06.03 Define future study sample sizes  ------------------------------------
n <- c(0, 100, seq(200, 2000, by = 200))
n.samples <- length(n)

### 06.04 Calculate EVSI for each parameter individually  ----------------------
# Initialize EVSI matrix for each parameter
evsi <- data.frame(N = n, matrix(0, nrow = n.samples, ncol = n.params))

# Name columns of EVSI matrix with parameter names
colnames(evsi)[-1] <- names.params[sel.params]

# Compute EVSI for all parameters separately
for (p in 1:n.params) { # p <- 1
  print(paste("Computing EVSI of parameter", names.params[sel.params[p]]))
  
  # Update loss based on Gaussian approximation for each sample of interest
  for (nSamp in 1:n.samples) { # nSamp <- 10
    Ltilde1 <- predict.ga(lmm1[[sel.params[p]]], n = n[nSamp], n0 = n0[p])
    Ltilde2 <- predict.ga(lmm2[[sel.params[p]]], n = n[nSamp], n0 = n0[p])
    
    # Combine losses into one matrix
    Ltilde <- cbind(Ltilde1, Ltilde2)
    
    # Apply EVSI equation
    evsi[nSamp, p + 1] <- mean(rowMaxs(Ltilde))
  }
}

### 06.05 Visualize EVSI results  ----------------------------------------------
# Create EVSI data frame for plotting in decreasing order of EVPPI
evsi.gg <- melt(evsi, 
                id.vars = "N", 
                variable.name = "Parameter", 
                value.name = "evsi")

evsi.gg$Parameter <- factor((evsi.gg$Parameter), 
                            levels = names.params[order(evppi.splines.gg$EVPPI, 
                                                        decreasing = TRUE)])

# Plot EVSI
ggplot(evsi.gg, aes(x = N, y = evsi)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Parameter) +
  ggtitle("Expected Value of Sample Information (EVSI)") +
  xlab("Sample size (n)") +
  ylab("$") +
  scale_x_continuous(breaks = number_ticks(5)) + 
  scale_y_continuous(breaks = number_ticks(6), labels = dollar) + 
  theme_bw(base_size = 14)

# Plot EVSI with EVPPI reference line
ggplot(evsi.gg, aes(x = N, y = evsi)) +
  geom_line(aes(linetype = "EVSI")) +
  geom_point() +
  facet_wrap(~ Parameter) +
  geom_hline(aes(yintercept = EVPPI, linetype = "EVPPI"), 
             data = evppi.splines.gg[sel.params, ]) +
  scale_linetype_manual(name = "", 
                        values = c("EVSI" = "solid", "EVPPI" = "dashed")) +
  xlab("Sample size (n)") +
  ylab("$") +
  scale_x_continuous(breaks = number_ticks(5)) + 
  scale_y_continuous(breaks = number_ticks(6), labels = dollar) + 
  theme_bw(base_size = 14)
