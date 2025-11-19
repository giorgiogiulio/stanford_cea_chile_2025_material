# *****************************************************************************
#
# Script: cSTM_sick-sicker_time_exercise_template.R
#
# Purpose: Implement and test the simulation-time-dependent Markov Sick-Sicker 
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
# This code implements a simulation-time-dependent Sick-Sicker cSTM model to 
# conduct a CEA of two strategies:
# - Standard of Care (SoC): best available care for the patients with the disease. 
#   This scenario reflects the natural history of the disease progression.
# - Strategy AB: This strategy combines treatment A and treatment B. The strategy 
#   treats both those Sick and Sicker. However, only for those Sick the treatment 
#   has an effect. For Sick individuals the disease progression is reduced, and 
#   individuals in the Sick state have an improved quality of life.
#
# *****************************************************************************

# ******************************************************************************
# 01 Exercise ----------------------------------------------------------
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
  #   Under the treatment strategy, individuals in the sick and sicker states are 
  #   treated until they recover (only if sick; individuals in the sicker state cannot 
  #   recover) or die. The cost of the treatment is additive to the baseline healthcare 
  #   costs of being sick or sicker. The treatment improves quality of life for those 
  #   individuals who are sick but has no impact on the quality of life of those who 
  #   are sicker. Unfortunately, it is not possible to reliably differentiate between 
  #   people in the sick and sicker states, so treatment cannot be targeted to only 
  #   those in the sick state. You are asked to evaluate the cost-effectiveness of 
  #   the treatment.
  #
  #   To model this disease, we will rely on a state-transition cohort model, called 
  #   the Sick-Sicker model, first described by Enns et al. The Sick-Sicker model 
  #   consists of four health states: Healthy (H), two disease states, Sick (S1) and 
  #   Sicker (S2), and Dead (D) (Figure 1). All individuals start in the Healthy state. 
  #   Over time, healthy individuals may develop the disease and can progress to S1. 
  #   Individuals in S1 can recover (return to state H), progress further to S2 or die. 
  #   Individuals in S2 cannot recover (i.e. cannot transition to either S1 or H). 
  #   Individuals in H have a baseline probability of death; individuals in S1 and S2 
  #   experience increased mortality compared to those in the H state, given in terms 
  #   of hazard ratios. These ratios are used to calculate the probabilities of dying 
  #   when in S1 and S2.
#
# ******************************************************************************
# 02 Table ----------------------------------------------------------
# ******************************************************************************
### 02.01 Model parameters  ----------------------------------------------------
#
# |           Parameter                |  R name              |   Value         |
# |:-----------------------------------|:---------------------|:---------------:|
# | Cycle length                       | `cycle_length`       | 1 year          |
# | Age at baseline                    | `n_age_init`         | 25 years old    |
# | Maximum age of follow-up           | `n_age_max`          | 100 years old   |
# | Names of health states             | `v_names_states`     | H, S1, S2, D.   |
# | Names of cycles (time horizon)     | `n_cycles`           | (n_age_max - n_age_init) / 
#                                                                  cycle_length |

### 02.02 Discount rates  --------------------------------------------------------
# | Annual discount rate (costs/QALYs) | `d_c` / `d_e`        | 3%              |
# | Annual transition probabilities conditional on surviving  | 
# | - Rate of becoming S1 when H       | `r_HS1`              | 0.15            |
# | - Rate of becoming H when S1       | `r_S1H`              | 0.5             |
# | - Rate of becoming S2 when S1      | `r_S1S2`             | 0.105           |

### 02.03 Annual mortality  -------------------------------------------
# | - All-cause mortality (H to D)     | `v_r_HD`             | HMD - info below|
# | - Hazard ratio of death in S1 vs H | `hr_S1`              | 3               |
# | - Hazard ratio of death in S2 vs H | `hr_S2`              | 10              |
# | - Hazard ratio of becoming Sicker  | `hr_S1S2_trtAB`      | 0.6             |
# |   when Sick under Strategy AB      |                      |                 |

### 02.04 Annual costs  --------------------------------------------------------
# | - Healthy individuals              | `c_H`                | $2,000          |
# | - Sick individuals in S1           | `c_S1`               | $4,000          |
# | - Sick individuals in S2           | `c_S2`               | $15,000         |
# | - Dead individuals                 | `c_D`                | $0              |
# | - Additional costs of sick                                                  |
# |   individuals treated with           `c_trtAB`              $25,000         |
# |   Strategy AB in S1 or S2                                                   |

### 02.05 Utility weights  ---------------------------------------------------
# | - Healthy individuals              | `u_H`                | 1.00            |
# | - Sick individuals in S1           | `u_S1`               | 0.75            |
# | - Sick individuals in S2           | `u_S2`               | 0.50            |
# | - Dead individuals                 | `u_D`                | 0.00            |
# | - Utility for individuals treated 
#     with Strategy AB in S1           | `u_trtAB`            | 0.95            |

# *Note:*
#   To calculate the probability of dying from S1 and S2, use the hazard ratios 
#   provided. Multiply the rate of dying from healthy by the appropriate hazard 
#   ratio, then convert the rate back to a probability. You can convert between 
#   rates and probabilities using:
#       r = -log(1 - p)
#       p = 1 - exp(-r * t)
#   The `darthtools` package includes helper functions `prob_to_rate()` and 
#   `rate_to_prob()` for this purpose.
#
# *HMD (Human Mortality Database):*  
#   The file `HMD_USA_Mx_2015.csv` includes age-specific mortality rates with 
#   columns: X, Year, Age, Female, Male, Total, and OpenInterval.  
#   Use the overall age-specific mortality rate from the `Total` column.

# ******************************************************************************
# 03 Setup ----------------------------------------------------------
# ******************************************************************************

### 03.01 Load packages and clear memory  --------------------------------------

rm(list = ls())    # clear memory (removes all the variables from the workspace)

# Install pacman if not present
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# Load pacman
library(pacman)

# Load (install if needed) CRAN packages
p_load(
  dplyr, tidyr, devtools, scales, ellipse, ggplot2, ggrepel, gridExtra,
  igraph, truncnorm, ggraph, patchwork, knitr, stringr, diagram, dampack
)

# Load (install if needed) GitHub packages
p_load_gh("DARTH-git/darthtools")

# ******************************************************************************
# 04 Model inputs ----------------------------------------------------------
# ******************************************************************************

### 04.01 General setup  ------------------------------------------------------
cycle_length <- 1   # cycle length equal to one year (use 1/12 for monthly)
n_age_init   <- 25  # age at baseline
n_age_max    <- 100 # maximum age of follow up
n_cycles     <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles

### 04.02 Age Labels  ----------------------------------------------------------
v_age_names  <- paste(rep(n_age_init:(n_age_max-1), each = 1/cycle_length), 
                      1:(1/cycle_length), 
                      sep = ".")

### 04.03 Health states  -------------------------------------------------------
# the 4 health states of the model:
v_names_states <- c("H",  # Healthy (H)
                    "S1", # Sick (S1)
                    "S2", # Sicker (S2)
                    "D")  # Dead (D)

n_states <- length(v_names_states)   # number of health states 

### 04.04 Discounting factors  -------------------------------------------------

# Comment out/ for the relevant country
# USA
d_c <- d_e <- 0.03  # annual discount rate for costs and QALY
# Canada
#d_c <- d_e <- 0.015 # annual discount rate for costs 

# NL
#d_c <- 0.03         # annual discount rate for costs 
#d_e <- 0.015        # annual discount rate for QALYs

### 04.05 Strategies  ------------------------------------------------------
v_names_str <- c("Standard of care", # store the strategy names
                 "Strategy AB") 
n_str       <- length(v_names_str)   # number of strategies

## Within-cycle correction (WCC) using Simpson's 1/3 rule 
v_wcc  <- gen_wcc(n_cycles = n_cycles, method = "Simpson1/3")

### 04.06 Transition rates  ----------------------------------------------------

### Transition rates (annual), and hazard ratios (HRs) 
r_HS1  <- 0.15  # constant annual rate of becoming Sick when Healthy
r_S1H  <- 0.5   # constant annual rate of becoming Healthy when Sick
r_S1S2 <- 0.105 # constant annual rate of becoming Sicker when Sick
hr_S1  <- 3     # hazard ratio of death in Sick vs Healthy 
hr_S2  <- 10    # hazard ratio of death in Sicker vs Healthy 

### 04.07 Treatment effect  ----------------------------------------------------

### Effectiveness of treatment AB 
hr_S1S2_trtAB <- 0.6  # hazard ratio of becoming Sicker when Sick under treatment AB

# ******************************************************************************
# 05 Load Age-dependent mortality rates data  ----------------------------------
# ******************************************************************************

## Age-dependent mortality rates 
lt_usa_2015 <- read.csv("data/HMD_USA_Mx_2015.csv")

### 05.01 Extract age-specific all-cause mortality rates  ----------------------
# Extract age-specific all-cause mortality for ages in model time horizon
v_r_mort_by_age <- lt_usa_2015 %>% 
  dplyr::filter(Age >= n_age_init & Age < n_age_max) %>%
  dplyr::select(Total) %>%
  as.matrix()

# ******************************************************************************
# 06 State rewards  ------------------------------------------------------------
# ******************************************************************************

### 06.01 Define state costs and utilities  ------------------------------------
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

### 06.02 Create vectors of state costs and utilities  --------------------------
### Discount weight for costs and effects 
v_dwc   <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
v_dwe   <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))

### 06.03 Process model inputs  ----------------------------------------------------
## Age-specific transition rates to the Dead state for all cycles 
v_r_HD_age  <- rep(v_r_mort_by_age, each = 1/cycle_length)
# Name age-specific mortality vector 
names(v_r_HD_age) <- v_age_names


### 06.04 Compute mortality rates  ----------------------------------------------
v_r_S1D_age <- v_r_HD_age * hr_S1 # Age-specific mortality rate in the Sick state 
v_r_S2D_age <- v_r_HD_age * hr_S2 # Age-specific mortality rate in the Sicker state 

### 06.05 Compute transition probabilities  ------------------------------------
# transform rates to probabilities adjusting by cycle length

# Constant annual probability of becoming Sick when Healthy conditional on surviving 
p_HS1       <- rate_to_prob(r = r_HS1,   t = cycle_length) 

# Constant annual probability of becoming Healthy when Sick conditional on surviving
p_S1H       <- rate_to_prob(r = r_S1H,   t = cycle_length) 

# Constant annual probability of becoming Sicker when Sick conditional on surviving
p_S1S2      <- rate_to_prob(r = r_S1S2,  t = cycle_length) 

# Age-specific mortality risk in the Healthy state 
v_p_HD_age  <- rate_to_prob(v_r_HD_age,  t = cycle_length) 

# Age-specific mortality risk in the Sick state
v_p_S1D_age <- rate_to_prob(v_r_S1D_age, t = cycle_length) 

# Age-specific mortality risk in the Sicker state
v_p_S2D_age <- rate_to_prob(v_r_S2D_age, t = cycle_length) 

### 06.06 Treatment effect on transition probabilities  --------------------------
## Annual transition probability of becoming Sicker when Sick for treatment AB 
# Apply hazard ratio to rate to obtain transition rate of becoming Sicker when Sick for treatment AB
r_S1S2_trtAB <- r_S1S2 * hr_S1S2_trtAB
# Transform rate to probability to become Sicker when Sick under treatment AB 
# adjusting by cycle length conditional on surviving
p_S1S2_trtAB <- rate_to_prob(r = r_S1S2_trtAB, t = cycle_length)


# ******************************************************************************
# 07 Exercise - Task 1 ------------------------------------------------------
# ******************************************************************************

# Build the Markov model in `R` for Standard of Care (SoC) and Strategy AB 
#  including age-dependent mortality and do the following:
# (1) Run the model
# (2) Plot the cohort trace from the model
# (3) Compute state rewards and expected outcomes (total utilities and costs)

# All starting healthy
v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
v_m_init

### 07.01 Initialize cohort traces  --------------------------------------------

### Initialize cohort trace under SoC 
m_M_SoC <- matrix(NA, 
              nrow = (n_cycles + 1), ncol = n_states, 
              dimnames = list(0:n_cycles, v_names_states))
# Store the initial state vector in the first row of the cohort trace
m_M_SoC[1, ] <- v_m_init


### 07.02 Initialize cohort trace for strategy AB  -----------------------------

# Structure and initial states are the same as for SoC
m_M_strAB <- m_M_SoC # Strategy AB

# ******************************************************************************
# 08 YOUR TURN: Create transition probability arrays  ---------------------------
# ******************************************************************************


### 08.01 Create transition probability arrays  --------------------------------
### 08.02 Create transition probability arrays for strategy SoC  ---------------

# # All transitions to a non-death state are assumed to be conditional on survival
 a_P_SoC <- array(0,
                  dim  = c(n_states, n_states, n_cycles),
                  dimnames = list(v_names_states,
                                  v_names_states,
                                  0:(n_cycles - 1)))


### 08.03 Fill in array From H -----------------

# your turn

### 08.04 Fill in array From S1 ----------------

# your turn

### 08.05 Fill in array From S2 ----------------

# your turn

### 08.06 Fill in array From D -----------------

# your turn

### 08.07 Initialize transition probability array for strategy AB  -------------

# your turn

### 08.08 Update only transition probabilities from S1 involving p_S1S2 --------

# your turn

### 08.09 CHECKS: ------------------------------------------------------------

# Check if transition probability arrays are valid
# Check that transition probabilities are [0, 1]

check_transition_probability(a_P_SoC,   verbose = TRUE)
check_transition_probability(a_P_strAB, verbose = TRUE)

# Check that all rows for each slice of the array sum to 1

check_sum_of_transition_array(a_P_SoC,   n_states = n_states, 
                              n_cycles = n_cycles, verbose = TRUE)

check_sum_of_transition_array(a_P_strAB, n_states = n_states, 
                              n_cycles = n_cycles, verbose = TRUE)

# ******************************************************************************
# 09 YOUR TURN: Implement the Markov model  ------------------------------------
# ******************************************************************************

### 09.01 Iterative solution of age-dependent cSTM  ----------------------------

# your turn

### 09.02 Store the cohort traces in a list   ----------------------------------

# your turn

# ******************************************************************************
# 10 YOUR TURN: Plot Outputs  --------------------------------------------------
# ******************************************************************************

### 10.01 Plot the cohort trace for strategies SoC and AB  ---------------------

# your turn
# ADD AN DESCRIPTION OF WHAT YOU SEE AND IF THAT MAKES SENSE TO YOU 

# ******************************************************************************
# 11 State Rewards   ----------------------------------------------------------
# ******************************************************************************

### 11.01 Scale by the cycle length  -------------------------------------------

### 11.02 Vector of state utilities under strategy SoC  ------------------------

# your turn

### 11.03 Vector of state costs under strategy SoC  -----------------------------

# your turn

### 11.06 Vector of state utilities under strategy AB  --------------------------

# your turn

### 11.07 Vector of state costs under strategy AB -------------------------------

# your turn


### 11.08 Store state rewards  -------------------------------------------------
# Store the vectors of state utilities for each strategy in a list 

# your turn

### 11.09  Store state costs  ----------------------------------------------------
# Store the vectors of state cost for each strategy in a list 

# your turn

### 11.10  Assign strategy names  ----------------------------------------------
# assign strategy names to matching items in the lists

# your turn


# ******************************************************************************
# 12 YOUR TURN: Compute expected outcomes   -----------------------------------
# ******************************************************************************


### 12.01 Calculate total utilities and costs  ---------------------------------
# Create empty vectors to store total utilities and costs 

# your turn

### 12.02 Loop through each strategy and calculate total utilities and costs  --

# your turn

#ADD AN INTERPRETATION OF YOUR RESULTS


# ******************************************************************************
# 12 Exercise - Task 2   -------------------------------------------------------
# ******************************************************************************

### 12.01 Estimate the cost-effectiveness of Strategy AB vs SoC.  --------------

### 12.02 Cost-effectiveness analysis (CEA) ------------------------------------ 

### 12.03 Incremental cost-effectiveness ratios (ICERs) ------------------------

# your turn

# ADD AN INTERPRETATION OF YOUR RESULTS


# ******************************************************************************
# 13 Exercise - Task 3   -------------------------------------------------------
# ******************************************************************************

### 13.01 Prepare results for CEA table  -----------------------------------------
#Create a well-formatted cost-effectiveness table with all results of interest and 
# plot the cost-effectiveness frontier.

# your turn

# ADD AN INTERPRETATION OF YOUR RESULTS

### 13.02 CEA frontier --------------------------------------------------------

# your turn

# ADD AN INTERPRETATION OF YOUR RESULTS

# *****************************************************************************
# 14 Acknowledgements  ---------------------------------------------------------
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
# End of cSTM_sick-sicker_time_exercise_template.R
# *****************************************************************************