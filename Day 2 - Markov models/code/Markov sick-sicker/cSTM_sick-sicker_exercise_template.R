# *****************************************************************************
#
# Script: cSTM_sick-sicker_exercise_template.R
#
# Purpose: Implement and test the time-independent Markov Sick-Sicker 
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
# 07 YOUR TURN: Initialize cohort traces  -------------------------------------
# ******************************************************************************

### 07.01 Initialize cohort trace for SoC  -------------------------------------

# your turn

### 07.02 Initialize cohort trace for strategy AB  -----------------------------
# Structure and initial states are the same as for SoC

# your turn

# ******************************************************************************
# 08 YOUR TURN: Create transition probability matrices  -----------------------
# ******************************************************************************

### 08.01 Create transition probability matrices for strategy SoC  -------------
### Initialize transition probability matrix for strategy SoC 
# All transitions to a non-death state are assumed to be conditional on survival 

# your turn

### 08.02 Fill in matrix  ------------------------------------------------------
# From H

# your turn

# From S1

# your turn

# From S2

# your turn

# From D

# your turn

### 08.03 Initialize transition probability matrix for strategy AB  ------------

# your turn

### 08.04 Update only transition probabilities from S1 involving p_S1S2  -------

# your turn

### 08.05 Check if transition probability matrices are valid  ------------------
### Check that transition probabilities are [0, 1] 

# your turn

### Check that all rows sum to 1 

# your turn

# ******************************************************************************
# 09 YOUR TURN: Run Markov model  ----------------------------------------------
# ******************************************************************************

### 09.01 Iterative solution of time-independent cSTM  -------------------------

# your turn

### 09.02 Store the cohort traces in a list  -----------------------------------

# your turn

# ******************************************************************************
# 10 YOUR TURN: Plot Outputs  --------------------------------------------------
# ******************************************************************************

### 10.01 Plot the cohort trace for strategies SoC and AB  ---------------------

# your turn

# ADD A DESCRIPTION OF WHAT YOU SEE AND IF THAT MAKES SENSE TO YOU 

# ******************************************************************************
# 11 YOUR TURN: State Rewards  -------------------------------------------------
# ******************************************************************************

### 11.01 Scale by the cycle length  -------------------------------------------

### 11.02 Vector of state utilities under strategy SoC  ------------------------

# your turn

### 11.03 Vector of state costs under strategy SoC  ----------------------------

# your turn

### 11.04 Vector of state utilities under strategy AB  -------------------------

# your turn

### 11.05 Vector of state costs under strategy AB  -----------------------------

# your turn

### 11.06 Store state rewards  -------------------------------------------------
# Store the vectors of state utilities for each strategy in a list 

# your turn

# Store the vectors of state cost for each strategy in a list 

# your turn

# assign strategy names to matching items in the lists

# your turn

# ******************************************************************************
# 12 YOUR TURN: Compute expected outcomes  ------------------------------------
# ******************************************************************************

### 12.01 Create empty vectors to store total utilities and costs  -------------

# your turn

### 12.02 Loop through each strategy and calculate total utilities and costs  --

# your turn

# ADD AN INTERPRETATION OF YOUR RESULTS

# ******************************************************************************
# 13 Exercise - Task 2  --------------------------------------------------------
# ******************************************************************************

### 13.01 Estimate the cost-effectiveness of Strategy AB vs SoC  ---------------

### 13.02 Cost-effectiveness analysis (CEA)  -----------------------------------

### 13.03 Incremental cost-effectiveness ratios (ICERs)  -----------------------

# your turn

# ADD AN INTERPRETATION OF YOUR RESULTS

# ******************************************************************************
# 14 Exercise - Task 3  --------------------------------------------------------
# ******************************************************************************

### 14.01 CEA table in proper format  ------------------------------------------
# Create a well-formatted cost-effectiveness table with all results of interest

# your turn

# ADD AN INTERPRETATION OF YOUR RESULTS

### 14.02 CEA frontier  --------------------------------------------------------
# Plot the cost-effectiveness frontier

# your turn

# ADD AN INTERPRETATION OF YOUR RESULTS

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
# End of cSTM_sick-sicker_exercise_template.R
# *****************************************************************************