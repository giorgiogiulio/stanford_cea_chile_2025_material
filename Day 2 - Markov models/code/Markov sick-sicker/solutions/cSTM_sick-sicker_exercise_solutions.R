# *****************************************************************************
#
# Script: cSTM_sick-sicker_exercise_solutions.R
#
# Purpose: SOLUTIONS - Implement and test the time-independent Markov Sick-Sicker 
#          cohort state-transition model for cost-effectiveness analysis (CEA).
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
#
# *****************************************************************************
#
# Notes:
#
# Please acknowledge our work. See details to cite below. 
#
# This code implements a time-independent Sick-Sicker cSTM model to 
# conduct a CEA of two strategies:
# - Standard of Care (SoC): best available care for the patients with the disease. 
#   This scenario reflects the natural history of the disease progression.
# - Strategy AB: This strategy combines treatment A and treatment B. The strategy 
#   treats both those Sick and Sicker. For Sick individuals the disease progression 
#   is reduced, and individuals in the Sick state have an improved quality of life.
#
# *****************************************************************************

# ******************************************************************************
# 01 Exercise ------------------------------------------------------------------
# ******************************************************************************

### 01.01 Instructions  --------------------------------------------------------

# Exercise: Construct a Markov Model of the Sick-Sicker Disease
#   In this exercise, we will model a hypothetical disease that affects individuals 
#   with an average age of 25 years and results in increased mortality, increased 
#   healthcare costs, and reduced quality of life. The disease has two levels; 
#   affected individuals initially become sick but can subsequently progress and 
#   become sicker. Two alternative strategies exist for this hypothetical disease: 
#   Standard of Care (SoC) and a treatment strategy - Strategy AB. 
#
#   Under the treatment strategy (Strategy AB), individuals in both the sick and 
#   sicker states receive treatment until they recover (only if sick; individuals 
#   in the sicker state cannot recover) or die. The cost of the treatment is 
#   additive to the baseline healthcare costs of being sick or sicker. The treatment 
#   provides two benefits: it improves quality of life for individuals who are sick 
#   and reduces the rate of disease progression from sick to sicker. However, it is 
#   not possible to reliably differentiate between people in the sick and sicker 
#   states, so treatment must be given to all individuals in both disease states and 
#   cannot be targeted to only those who would benefit most. You are asked to 
#   evaluate the cost-effectiveness of Strategy AB compared to Standard of Care.
#
#   To model this disease, we will rely on a time-independent state-transition 
#   cohort model, called the Sick-Sicker model, first described by Enns et al. 
#   The Sick-Sicker model consists of four health states: Healthy (H), two disease 
#   states, Sick (S1) and Sicker (S2), and Dead (D). All individuals start in the 
#   Healthy state at age 25 and are followed until age 100. Over time, healthy 
#   individuals may develop the disease and progress to S1. Individuals in S1 can 
#   recover (return to state H), progress further to S2, or die. Individuals in S2 
#   cannot recover (i.e., cannot transition to either S1 or H) and can only remain 
#   in S2 or die. Individuals in H have a baseline probability of death (all-cause 
#   mortality); individuals in S1 and S2 experience increased mortality compared to 
#   those in the H state, given in terms of hazard ratios. These ratios are used to 
#   calculate the annual probabilities of dying when in S1 and S2.
#
#   Under Strategy AB, the hazard ratio for disease progression from S1 to S2 is 
#   reduced by 40% (hr_S1S2_trtAB = 0.6), reflecting the treatment's effectiveness 
#   in slowing disease progression. Additionally, individuals in S1 who receive 
#   treatment AB experience improved quality of life (utility = 0.95) compared to 
#   those receiving standard of care (utility = 0.75), while treatment has no effect 
#   on the quality of life of individuals in S2.
#
#   The model uses annual cycles with appropriate discounting (3% for both costs 
#   and QALYs) and applies within-cycle correction using Simpson's 1/3 rule. The 
#   analysis evaluates total costs and quality-adjusted life years (QALYs) for each 
#   strategy and calculates incremental cost-effectiveness ratios (ICERs) to 
#   determine the value of Strategy AB compared to Standard of Care.
#

# ******************************************************************************
# 02 Table ---------------------------------------------------------------------
# ******************************************************************************
### 02.01 Model parameters  ----------------------------------------------------
#
# |           Parameter                |  R name              |   Value         |
# |:-----------------------------------|:---------------------|:---------------:|
# | Cycle length                       | `cycle_length`       | 1 year          |
# | Age at baseline                    | `n_age_init`         | 25 years old    |
# | Maximum age of follow-up           | `n_age_max`          | 100 years old   |
# | Names of health states             | `v_names_states`     | H, S1, S2, D    |
# | Time horizon (number of cycles)    | `n_cycles`           | (n_age_max - n_age_init) / 
#                                                                  cycle_length |

### 02.02 Discount rates  ------------------------------------------------------
# | Annual discount rate (costs/QALYs) | `d_c` / `d_e`        | 3%              |

### 02.03 Transition rates  ----------------------------------------------------
# | Annual transition rates conditional on surviving:                           | 
# | - Rate of becoming S1 when H       | `r_HS1`              | 0.15            |
# | - Rate of becoming H when S1       | `r_S1H`              | 0.5             |
# | - Rate of becoming S2 when S1      | `r_S1S2`             | 0.105           |

### 02.04 Mortality rates  -----------------------------------------------------
# | - All-cause mortality rate (H to D)| `r_HD`               | 0.002           |
# | - Hazard ratio of death in S1 vs H | `hr_S1`              | 3               |
# | - Hazard ratio of death in S2 vs H | `hr_S2`              | 10              |

### 02.05 Treatment effectiveness  ---------------------------------------------
# | - Hazard ratio of becoming Sicker  | `hr_S1S2_trtAB`      | 0.6             |
# |   when Sick under Strategy AB      |                      |                 |

### 02.06 Annual costs  --------------------------------------------------------
# | - Healthy individuals              | `c_H`                | $2,000          |
# | - Sick individuals in S1           | `c_S1`               | $4,000          |
# | - Sick individuals in S2           | `c_S2`               | $15,000         |
# | - Dead individuals                 | `c_D`                | $0              |
# | - Additional costs of treatment AB | `c_trtAB`            | $25,000         |
# |   for individuals in S1 or S2      |                      |                 |

### 02.07 Utility weights  -----------------------------------------------------
# | - Healthy individuals              | `u_H`                | 1.00            |
# | - Sick individuals in S1           | `u_S1`               | 0.75            |
# | - Sick individuals in S2           | `u_S2`               | 0.50            |
# | - Dead individuals                 | `u_D`                | 0.00            |
# | - Utility for individuals in S1    | `u_trtAB`            | 0.95            |
# |   treated with Strategy AB         |                      |                 |

# *Note:*
#   To calculate the probability of dying from S1 and S2, use the hazard ratios 
#   provided. Multiply the rate of dying from healthy by the appropriate hazard 
#   ratio, then convert the rate back to a probability using the formula:
#       p = 1 - exp(-r * t)
#   The `darthtools` package includes the helper function `rate_to_prob()` 
#   for this purpose.

# ******************************************************************************
# 03 Setup ---------------------------------------------------------------------
# ******************************************************************************

### 03.01 Load packages and clear memory  --------------------------------------

rm(list = ls())    # clear memory (removes all the variables from the workspace)

# Install pacman if not present
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# Load pacman
library(pacman)

# Load (install if needed) CRAN packages
p_load(
  dplyr, tidyr, devtools, scales, ellipse, ggplot2, lazyeval, igraph, 
  truncnorm, ggraph, reshape2, knitr, stringr, diagram, dampack
)

# Load (install if needed) GitHub packages
p_load_gh("DARTH-git/darthtools")

# ******************************************************************************
# 04 Model inputs --------------------------------------------------------------
# ******************************************************************************

### 04.01 General setup  -------------------------------------------------------
cycle_length <- 1   # cycle length equal to one year (use 1/12 for monthly)
n_age_init   <- 25  # age at baseline
n_age_max    <- 100 # maximum age of follow up
n_cycles     <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles

### 04.02 Health states  -------------------------------------------------------
# the 4 health states of the model:
v_names_states <- c("H",  # Healthy (H)
                    "S1", # Sick (S1)
                    "S2", # Sicker (S2)
                    "D")  # Dead (D)
# NOTE: For our parameter values of costs and utilities we use
# just letters for the health states 
# Healthy (H), Sick (S1), Sicker (S2), Dead (D)                                          

n_states <- length(v_names_states)   # number of health states 

### 04.03 Discounting factors  -------------------------------------------------
d_c <- 0.03        # annual discount rate for costs 
d_e <- 0.03        # annual discount rate for QALYs

### 04.04 Strategies  ----------------------------------------------------------
v_names_str <- c("Standard of care",  # store the strategy names
                 "Strategy AB") 
n_str       <- length(v_names_str)    # number of strategies

### 04.05 Within-cycle correction (WCC) using Simpson's 1/3 rule  --------------
v_wcc <- gen_wcc(n_cycles = n_cycles,  method = "Simpson1/3")

### 04.06 Transition rates (annual), and hazard ratios (HRs)  ------------------
r_HD    <- 0.002 # constant annual rate of dying when Healthy (all-cause mortality)
r_HS1   <- 0.15  # constant annual rate of becoming Sick when Healthy
r_S1H   <- 0.5   # constant annual rate of becoming Healthy when Sick
r_S1S2  <- 0.105 # constant annual rate of becoming Sicker when Sick
hr_S1   <- 3     # hazard ratio of death in Sick vs Healthy 
hr_S2   <- 10    # hazard ratio of death in Sicker vs Healthy 

### 04.07 Effectiveness of treatment AB  ---------------------------------------
hr_S1S2_trtAB <- 0.6  # hazard ratio of becoming Sicker when Sick under treatment AB

### 04.08 State rewards  -------------------------------------------------------
#### Costs 
c_H     <- 2000  # annual cost of being Healthy
c_S1    <- 4000  # annual cost of being Sick
c_S2    <- 15000 # annual cost of being Sicker
c_D     <- 0     # annual cost of being dead
c_trtAB <- 25000 # annual cost of receiving treatment AB
#### Utilities 
u_H     <- 1     # annual utility of being Healthy
u_S1    <- 0.75  # annual utility of being Sick
u_S2    <- 0.5   # annual utility of being Sicker
u_D     <- 0     # annual utility of being dead
u_trtAB <- 0.95  # annual utility when receiving treatment AB

### 04.09 Discount weight for costs and effects  -------------------------------
v_dwc   <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
v_dwe   <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))

### 04.10 Process model inputs  ------------------------------------------------
## Cycle-specific transition probabilities to the Dead state 
# compute mortality rates
r_S1D   <- r_HD * hr_S1 # annual mortality rate in the Sick state
r_S2D   <- r_HD * hr_S2 # annual mortality rate in the Sicker state

# transform rates to probabilities 
p_HS1   <- rate_to_prob(r = r_HS1,  t = cycle_length) # constant annual probability of becoming Sick when Healthy conditional on surviving 
p_S1H   <- rate_to_prob(r = r_S1H,  t = cycle_length) # constant annual probability of becoming Healthy when Sick conditional on surviving
p_S1S2  <- rate_to_prob(r = r_S1S2, t = cycle_length) # constant annual probability of becoming Sicker when Sick conditional on surviving
p_HD    <- rate_to_prob(r = r_HD,   t = cycle_length) # annual mortality risk in the Healthy state
p_S1D   <- rate_to_prob(r = r_S1D,  t = cycle_length) # annual mortality risk in the Sick state
p_S2D   <- rate_to_prob(r = r_S2D,  t = cycle_length) # annual mortality risk in the Sicker state

## Annual transition probability of becoming Sicker when Sick for treatment AB 
# Apply hazard ratio to rate to obtain transition rate of becoming Sicker when 
# Sick for treatment AB
r_S1S2_trtAB <- r_S1S2 * hr_S1S2_trtAB
# Transform rate to probability to become Sicker when Sick under treatment AB conditional on surviving
p_S1S2_trtAB <- rate_to_prob(r = r_S1S2_trtAB, t = cycle_length) 

# ******************************************************************************
# 05 Construct state-transition models -----------------------------------------
# ******************************************************************************

### 05.01 Create state transition diagram  -------------------------------------
m_P_diag <- matrix(0, nrow = n_states, ncol = n_states, 
                   dimnames = list(v_names_states, v_names_states))
m_P_diag["H" , "S1"] = "" 
m_P_diag["H" , "D" ] = "" 
m_P_diag["H" , "H" ] = "" 
m_P_diag["S1", "H" ] = "" 
m_P_diag["S1", "S2"] = "" 
m_P_diag["S1", "D" ] = "" 
m_P_diag["S1", "S1"] = "" 
m_P_diag["S2", "D" ] = "" 
m_P_diag["S2", "S2"] = "" 
m_P_diag["D", "D"  ] = "" 
layout.fig <- c(3, 1)

plotmat(t(m_P_diag), t(layout.fig), self.cex = 0.5, curve = 0, arr.pos = 0.7,  
        latex = T, arr.type = "curved", relsize = 0.9, box.prop = 0.8, 
        cex = 0.8, box.cex = 0.9, lwd = 1)

# ******************************************************************************
# 06 Exercise - Task 1 ---------------------------------------------------------
# ******************************************************************************

# Build the Markov model in `R` for Standard of Care (SoC) and Strategy AB 
# and do the following:
# (1) Initialize the cohort trace
# (2) Create transition probability matrices
# (3) Run the model
# (4) Plot the cohort trace from the model
# (5) Compute state rewards and expected outcomes (total utilities and costs)

### 06.01 Initial state vector  ------------------------------------------------
# All starting healthy
v_m_init <- c(Healthy = 1, Sick = 0, Sicker = 0, Dead = 0) # initial state vector
v_m_init

# ******************************************************************************
# 07 Initialize cohort traces  -------------------------------------------------
# ******************************************************************************

### 07.01 Initialize cohort trace for SoC  -------------------------------------
m_M <- matrix(NA, 
              nrow = (n_cycles + 1), ncol = n_states, 
              dimnames = list(0:n_cycles, v_names_states))
# Store the initial state vector in the first row of the cohort trace
m_M[1, ] <- v_m_init

### 07.02 Initialize cohort trace for strategy AB  -----------------------------
# Structure and initial states are the same as for SoC
m_M_strAB <- m_M # Strategy AB

# ******************************************************************************
# 08 Create transition probability matrices  -----------------------------------
# ******************************************************************************

### 08.01 Create transition probability matrices for strategy SoC  -------------
### Initialize transition probability matrix for strategy SoC 
# All transitions to a non-death state are assumed to be conditional on survival 

v_names_states <- c("Healthy",  # Healthy (H)
                    "Sick", # Sick (S1)
                    "Sicker", # Sicker (S2)
                    "Dead")  # Dead (D)

m_P <- matrix(0, 
              nrow = n_states, ncol = n_states, 
              dimnames = list(v_names_states, 
                              v_names_states)) # define row and column names

### 08.02 Fill in matrix  ------------------------------------------------------
# From H
m_P["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS1)
m_P["Healthy", "Sick"]    <- (1 - p_HD) *      p_HS1 
m_P["Healthy", "Dead"]    <-      p_HD

# From S1
m_P["Sick", "Healthy"]   <- (1 - p_S1D) *       p_S1H
m_P["Sick", "Sick"]      <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2))
m_P["Sick", "Sicker"]    <- (1 - p_S1D) *               p_S1S2
m_P["Sick", "Dead"]      <-      p_S1D

# From S2
m_P["Sicker", "Sicker"]  <- 1 - p_S2D
m_P["Sicker", "Dead"]    <-     p_S2D

# From D
m_P["Dead", "Dead"]      <- 1

### 08.03 Initialize transition probability matrix for strategy AB  ------------
m_P_strAB <- m_P

### 08.04 Update only transition probabilities from S1 involving p_S1S2  -------
# Update only transition probabilities from S1 involving p_S1S2
m_P_strAB["Sick", "Sick"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2_trtAB))
m_P_strAB["Sick", "Sicker"] <- (1 - p_S1D) *               p_S1S2_trtAB

### 08.05 Check if transition probability matrices are valid  ------------------
### Check that transition probabilities are [0, 1] 
check_transition_probability(m_P,       
                             verbose = TRUE) # m_P >= 0 && m_P <= 1
check_transition_probability(m_P_strAB, 
                             verbose = TRUE) # m_P_strAB >= 0 && m_P_strAB <= 1

### Check that all rows sum to 1 
check_sum_of_transition_array(m_P,
                              n_states = n_states, n_cycles = n_cycles, 
                              verbose = TRUE) # rowSums(m_P) == 1
check_sum_of_transition_array(m_P_strAB, 
                              n_states = n_states, n_cycles = n_cycles, 
                              verbose = TRUE) # rowSums(m_P_strAB) == 1

# ******************************************************************************
# 09 Run Markov model  ---------------------------------------------------------
# ******************************************************************************

### 09.01 Iterative solution of time-independent cSTM  -------------------------
# Iterative solution of time-independent cSTM
for (t in 1:n_cycles) {
  # For SoC
  m_M[t + 1, ] <- m_M[t, ] %*% m_P
  # For strategy AB
  m_M_strAB[t + 1, ] <- m_M_strAB[t, ] %*% m_P_strAB
}

### 09.02 Store the cohort traces in a list  -----------------------------------
## Store the cohort traces in a list 
l_m_M <- list(m_M,
              m_M_strAB)
names(l_m_M) <- v_names_str

# ******************************************************************************
# 10 Plot Outputs  -------------------------------------------------------------
# ******************************************************************************

### 10.01 Plot the cohort trace for strategies SoC and AB  ---------------------


#### 10.01.01 Plot SoC

# Convert to data frame
df_M <- data.frame(Cycle = 0:n_cycles, m_M, check.names = F)

# Convert to long format
df_M_long <- tidyr::gather(df_M, key = `Health State`, value, 
                           2:ncol(df_M))

# Assign factor levels to health states
df_M_long$`Health State` <- factor(
  df_M_long$`Health State`,
  levels = c("H", "S1", "S2", "D"),              # existing values
  labels = c("Healthy", "Sick", "Sicker", "Dead") # new labels
)

#Plot 
ggplot(df_M_long, aes(x = Cycle, y = value, color = `Health State`, 
                  linetype = `Health State`)) + 
  geom_line(size = 1) + 
  xlab("Cycle") + 
  ylab("Proportion of the cohort") + 
  scale_x_continuous(breaks = number_ticks(8)) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        legend.background = element_rect(fill = NA))


#### 10.01.02 Plot Strategy AB

# Convert to data frame

df_M_strAB <- data.frame(Cycle = 0:n_cycles, m_M_strAB, check.names = F)

# Convert to long format

df_M_strAB_long <- tidyr::gather(df_M_strAB, key = `Health State`, value, 
                               2:ncol(df_M_strAB))

# Assign factor levels to health states

df_M_strAB_long$`Health State` <- factor(
  df_M_strAB_long$`Health State`,
  levels = c("H", "S1", "S2", "D"),              # existing values
  labels = c("Healthy", "Sick", "Sicker", "Dead") # new labels
)

#Plot

ggplot(df_M_strAB_long, aes(x = Cycle, y = value, color = `Health State`, 
                          linetype = `Health State`)) + 
  geom_line(size = 1) + 
  xlab("Cycle") + 
  ylab("Proportion of the cohort") + 
  scale_x_continuous(breaks = number_ticks(8)) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        legend.background = element_rect(fill = NA))



# We combine both plots into one for easier comparison

# Combine data frames
df_M_long$strategy <- "Standard of Care"

df_M_strAB_long$strategy <- "Strategy AB"

df_strategies <- rbind(df_M_long, df_M_strAB_long)

# Plot combined

ggplot(df_strategies, aes(x = Cycle, y = value, color = `Health State`, 
                          linetype = `Health State`)) + 
  geom_line(size = 1) + 
  xlab("Cycle") + 
  ylab("Proportion of the cohort") + 
  scale_x_continuous(breaks = number_ticks(8)) + 
  facet_wrap(~ strategy) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        legend.background = element_rect(fill = NA))

# INTERPRETATION: 

# The cohort trace plots illustrate the distribution of the cohort across
# health states over time for both strategies.
#
# Under Standard of Care:

# - The proportion of individuals in the Healthy state decreases over time
#   as individuals transition to the Sick and Sicker states or die.
# - The Sick state initially increases as individuals become sick, but then
#   decreases as they either recover, progress to Sicker, or die.
# - The Sicker state increases over time as individuals progress from Sick
#   to Sicker, reflecting the disease progression.
# - The Dead state accumulates individuals over time as mortality occurs
#   from all health states.
#
# Under Strategy AB:

# - The proportion of individuals in the Healthy state remains higher
#   compared to Standard of Care, indicating that the treatment helps
#   maintain health.
# - The Sick state shows a slower increase and a lower peak compared to
#   Standard of Care, reflecting the treatment's effectiveness in improving
#   quality of life and reducing disease progression.
# - The Sicker state increases at a slower rate compared to Standard of Care,
#   indicating that the treatment effectively reduces progression from Sick
#   to Sicker.
# - The Dead state still accumulates individuals over time, but the overall
#   mortality may be lower due to the treatment's benefits.
# Overall, Strategy AB appears to improve health outcomes by maintaining
# a higher proportion of individuals in better health states and slowing disease
# progression compared to Standard of Care.


# ******************************************************************************
# 11 State Rewards  ------------------------------------------------------------
# ******************************************************************************

### 11.01 Scale by the cycle length  -------------------------------------------

### 11.02 Vector of state utilities under strategy SoC  ------------------------
# Vector of state utilities under strategy SoC
v_u_SoC    <- c(Healthy = u_H, 
                Sick    = u_S1, 
                Sicker  = u_S2, 
                Dead    = u_D) * cycle_length

### 11.03 Vector of state costs under strategy SoC  ----------------------------
# Vector of state costs under strategy SoC
v_c_SoC    <- c(Healthy = c_H, 
                Sick    = c_S1,
                Sicker  = c_S2, 
                Dead    = c_D) * cycle_length

### 11.04 Vector of state utilities under strategy AB  -------------------------
# Vector of state utilities under strategy AB
v_u_strAB  <- c(Healthy = u_H, 
                Sick    = u_trtAB, 
                Sicker  = u_S2, 
                Dead    = u_D) * cycle_length

### 11.05 Vector of state costs under strategy AB  -----------------------------
# Vector of state costs under strategy AB
v_c_strAB  <- c(Healthy = c_H, 
                Sick    = c_S1 + c_trtAB, 
                Sicker  = c_S2 + c_trtAB, 
                Dead    = c_D) * cycle_length

### 11.06 Store state rewards  -------------------------------------------------
## Store state rewards 
# Store the vectors of state utilities for each strategy in a list 
l_u   <- list(SQ = v_u_SoC,
              AB = v_u_strAB)

# Store the vectors of state cost for each strategy in a list 
l_c   <- list(SQ = v_c_SoC,
              AB = v_c_strAB)

# assign strategy names to matching items in the lists
names(l_u) <- names(l_c) <- v_names_str

# ******************************************************************************
# 12 Compute expected outcomes  ------------------------------------------------
# ******************************************************************************

### 12.01 Create empty vectors to store total utilities and costs  -------------
# Create empty vectors to store total utilities and costs 
v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str

### 12.02 Loop through each strategy and calculate total utilities and costs  --
# Loop through each strategy and calculate total utilities and costs 
for (i in 1:n_str) {
  v_u_str <- l_u[[i]] # select the vector of state utilities for the i-th strategy
  v_c_str <- l_c[[i]] # select the vector of state costs for the i-th strategy
  
  # Expected QALYs and costs per cycle 
  # Vector of QALYs and Costs
  # Apply state rewards 
  v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
  v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
  
  # Discounted total expected QALYs and Costs per strategy and apply 
  # within-cycle correction if applicable
  # QALYs
  v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
  # Costs
  v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
}

# Display results
v_tot_qaly
v_tot_cost

# INTERPRETATION OF RESULTS:
# Strategy AB results in:
# - Higher total QALYs compared to Standard of Care, due to:
#   * Improved quality of life for individuals in the Sick state (utility 0.95 vs 0.75)
#   * Reduced disease progression, meaning fewer individuals reach the Sicker state
#     with its lower utility (0.50)
# - Higher total costs, due to the additional treatment cost ($25,000 per year)
#   for all individuals in S1 and S2 states
# 
# The key question is whether the additional QALYs gained justify the additional
# costs incurred, which will be addressed through the cost-effectiveness analysis.

# ******************************************************************************
# 13 Exercise - Task 2  --------------------------------------------------------
# ******************************************************************************

### 13.01 Estimate the cost-effectiveness of Strategy AB vs SoC  ---------------

### 13.02 Cost-effectiveness analysis (CEA)  -----------------------------------

### 13.03 Incremental cost-effectiveness ratios (ICERs)  -----------------------
## Incremental cost-effectiveness ratios (ICERs) 
df_cea <- calculate_icers(cost       = v_tot_cost, 
                          effect     = v_tot_qaly,
                          strategies = v_names_str)
df_cea

# INTERPRETATION OF RESULTS:
# The ICER represents the additional cost per additional QALY gained with Strategy AB
# compared to Standard of Care. 
# 
# Decision-makers would compare this ICER to a willingness-to-pay threshold:
# - If the ICER is below the threshold, Strategy AB would be considered cost-effective
# - Common thresholds range from $50,000 to $150,000 per QALY in the US
# - Other countries use different thresholds (e.g., £20,000-30,000 per QALY in UK)
#
# The ICER tells us how much society would need to pay for each additional year
# of perfect health gained by implementing Strategy AB instead of Standard of Care.

# ******************************************************************************
# 14 Exercise - Task 3  --------------------------------------------------------
# ******************************************************************************

### 14.01 CEA table in proper format  ------------------------------------------
## CEA table in proper format 
table_cea <- format_table_cea(df_cea) 
table_cea

# INTERPRETATION:
# This formatted table presents the cost-effectiveness analysis results in a 
# standard format suitable for publication or reporting. It shows:
# - Total costs and QALYs for each strategy
# - Incremental costs and QALYs (the differences between strategies)
# - The ICER (incremental cost-effectiveness ratio)
# - Status: whether each strategy is on the cost-effectiveness frontier
#
# The table makes it easy to see both the absolute values for each strategy
# and the incremental differences that drive the cost-effectiveness decision.

### 14.02 CEA frontier  --------------------------------------------------------
## CEA frontier 
plot(df_cea, label = "all", txtsize = 14) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.85, 0.3))

# INTERPRETATION:
# The cost-effectiveness frontier plot visualizes the trade-off between costs
# and health outcomes (QALYs):
# - Each point represents a strategy
# - The x-axis shows total QALYs, the y-axis shows total costs
# - The slope of the line connecting strategies represents the ICER
# - Strategies on the frontier are potentially cost-effective options
# - Strategies below/right of the frontier are "dominated" (worse outcomes for more cost)
#
# This visual representation helps decision-makers understand:
# - The magnitude of additional costs required for additional health benefits
# - Whether any strategies are clearly dominated
# - The overall efficiency of each strategy in converting resources to health outcomes

# *****************************************************************************
# 15 Acknowledgements  ---------------------------------------------------------
# *****************************************************************************

# We kindly request you to add the following Acknowledgement paragraph to your
# further work where DARTH code formed the basis. You may also include additional
# sources of reference to acknowledge other contributors whose code you have used.
#
# For this work, we made use of the template developed by the
# Decision Analysis in R for Technologies in Health (DARTH) workgroup:
# http://darthworkgroup.com
#
# The notation of our code is based on the following framework and coding convention:
# Alarid-Escudero, F., Krijkamp, E., Pechlivanoglou, P. et al.
# "A Need for Change! A Coding Framework for Improving Transparency in Decision Modeling."
# PharmacoEconomics 37, 1329–1339 (2019).
# https://doi.org/10.1007/s40273-019-00837-x
#
# Other work from DARTH can be found at:
# http://darthworkgroup.com/publications/
#
# Copyright for Assignment Work
#
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS.
# All rights reserved in Canada, the United States, and worldwide.
# Copyright, trademarks, trade names, and any and all associated intellectual property
# are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the collaborating institutions.
#
# These materials may be used, reproduced, modified, distributed, and adapted
# with proper attribution.
# End of cSTM_sick-sicker_exercise_solutions.R
# ***