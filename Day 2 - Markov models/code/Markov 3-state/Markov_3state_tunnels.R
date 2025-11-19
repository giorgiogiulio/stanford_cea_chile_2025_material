#*******************************************************************************
# Developed by the Decision Analysis in R for Technologies in Health (DARTH) 
# workgroup:
#   
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (2)	
# M.G. Myriam Hunink, MD, PhD (3,4)
# Hawre J. Jalal, MD, PhD (5) 
# Eline M. Krijkamp, MSc (3)	
# Petros Pechlivanoglou, PhD (6)

# Affiliations: 		
# 1 Department of Health Policy, School of Medicine, and Stanford Health Policy, 
#   and Freeman-Spogli Institute for International Studies, Stanford University, 
#   Stanford, California, USA
# 2 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 3 Erasmus MC, Rotterdam, The Netherlands
# 4 Harvard T.H. Chan School of Public Health, Boston, USA
# 5 University of Ottawa, Ottawa, Ontario, Canada
# 6 The Hospital for Sick Children, Toronto and University of Toronto, Toronto, 
#   ON, Canada

# Please cite our publications when using this code:
# - Alarid-Escudero, F, Krijkamp, E, Enns, EA, Yang, A, Hunink, MGGM, 
#   Pechlivanoglou, P, & Jalal, H. (2023). An Introductory Tutorial on Cohort 
#   State-Transition Models in R Using a Cost-Effectiveness Analysis Example. 
#   Medical Decision Making, 43(1):3–20. 
#   https://doi.org/10.1177/0272989X221103163
#
# - Alarid-Escudero, F, Krijkamp, E, Enns, EA, Yang, A, Hunink, MGGM, 
#   Pechlivanoglou, P, & Jalal, H. (2023). A Tutorial on Time-Dependent Cohort 
#   State-Transition Models in R using a Cost-Effectiveness Analysis Example. 
#   Medical Decision Making, 43(1):21–41. 
#   https://doi.org/10.1177/0272989X221121747
#   
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, Hunink MG. 
#   (2017). An Overview of R in Health Decision Sciences. Med Decis Making, 
#   37(3):735-746. https://journals.sagepub.com/doi/abs/10.1177/0272989X16686559
#  
# - Krijkamp EM, Alarid-Escudero F, Enns E, Pechlivanoglou P, Hunink MM, Jalal H. 
#   (2020).A Multidimensional Array Representation of State-Transition Model 
#   Dynamics. Med Decis Making, 40(2):242-248 
#   https://doi.org/10.1177/0272989X19893973
# 
#*******************************************************************************

#******************************************************************************
### MARKOV MODEL - THREE STATE WITH TUNNELS ####
#******************************************************************************

# Clear memory (removes all the variables from the workspace)
rm(list = ls()) 

#******************************************************************************
#### 01 Load packages ####
#******************************************************************************

# No packages required

#******************************************************************************
##### 02 Load functions ####
#******************************************************************************

# No functions required

#******************************************************************************
##### 03 Input model parameters ####
#******************************************************************************

# Strategy names
v_names_str <- c("Base Case")  

# Number of strategies
n_str <- length(v_names_str)

# Markov model parameters
v_n      <- c("Healthy", "Sick", "Dead") # State names
n_states <- length(v_n)                  # Number of states
n_age_init <- 25  # age at baseline
n_age_max  <- 100 # maximum age of follow up
n_t        <- n_age_max - n_age_init # Number of cycles

# Tunnels
n_tunnel_size <- n_t

# Vector with cycles for tunnels
v_cycles_tunnel <- 1:n_tunnel_size

# Sick state
v_sick_tunnels <- paste("Sick_", seq(1, n_tunnel_size), "Yr", sep = "")

# Create variables for time-dependent model
v_n_tunnels       <- c("Healthy", v_sick_tunnels, "Dead")  # State names
n_states_tunnels  <- length(v_n_tunnels)                   # Number of states

v_p_HD <- seq(0.003, 0.01, length.out = n_t)  # Probability to die when sick (age-dependent)
p_HS <- 0.20                                  # Probability to become sick when healthy

#* Weibull parameters for state-residence-dependent transition probability of 
#* dying when Sick
r_SD_scale <- 0.08 # scale
r_SD_shape <- 1.1  # shape

## State-residence-dependent transition rate of dying when Sick ----
#* Weibull transition rate
v_r_SD_tunnels <- (v_cycles_tunnel*r_SD_scale)^r_SD_shape - 
  ((v_cycles_tunnel - 1)*r_SD_scale)^r_SD_shape
#* Weibull transition probability 
v_p_SD_tunnels <- 1 - exp(-v_r_SD_tunnels)

# Costs and utilities  
c_H  <- 400                                 # Cost of remaining one cycle healthy
c_S  <- 1000                                # Cost of remaining one cycle sick
c_D  <- 0                                   # Cost of remaining one cycle dead
u_H  <- 1                                   # Utility when healthy 
u_S  <- 0.5                                 # Utility when sick
u_D  <- 0                                   # Utility when dead
d_e <- d_c <- 0.03                          # Equal discount of costs and QALYs by 3%
v_costs     <- c(c_H, c_S, c_D)             # All costs
v_utilities <- c(u_H, u_S, u_D)             # All utilities

# Calculate discount weights for costs for each cycle based on discount rate d_c
v_dwc <- 1 / (1 + d_e) ^ (0:n_t) 

# Calculate discount weights for effectiveness for each cycle based on discount rate d_e
v_dwe <- 1 / (1 + d_c) ^ (0:n_t) 

#******************************************************************************
#### 04 Define and initialize matrices and vectors ####
#******************************************************************************
#******************************************************************************
#### 04.1 Cohort trace ####
#******************************************************************************

# Create Markov trace (n_t + 1 because R doesn't understand Cycle 0)
m_M <- matrix(NA, 
              nrow     = n_t + 1,  
              ncol     = n_states_tunnels,                  
              dimnames = list(0:n_t, v_n_tunnels)) 

# The cohort starts as healthy
# Initialize first cycle of Markov trace accounting for the tunnels
m_M[1, ] <- c(1, rep(0, n_tunnel_size), 0) 

#******************************************************************************
#### 04.2 Transition probability array ####
#******************************************************************************

# Create the transition probability array
a_P <- array(0,
             dim      = c(n_states_tunnels, n_states_tunnels, n_t),               
             dimnames = list(v_n_tunnels, v_n_tunnels, 0:(n_t - 1)))  

# Fill in the transition probability array:
# From Healthy
a_P["Healthy", "Healthy", ]  <- 1 - v_p_HD - p_HS
a_P["Healthy", "Sick_1Yr", ] <- p_HS
a_P["Healthy", "Dead", ]     <- v_p_HD

# From Sick
# From Sick to Sick in each step of the tunnel minus the last one
for (i in 1:(n_tunnel_size - 1)) { #i=1
  a_P[v_sick_tunnels[i], v_sick_tunnels[i + 1], ] <- 1 - v_p_SD_tunnels[i]
  a_P[v_sick_tunnels[i], "Dead", ] <- v_p_SD_tunnels[i]
}

# From Sick_60Yr to Sick Sick_60Yr
a_P[v_sick_tunnels[n_tunnel_size], 
    v_sick_tunnels[n_tunnel_size], ] <- 1 - v_p_SD_tunnels[n_tunnel_size]

# From Sick_60Yr to Death
a_P[v_sick_tunnels[n_tunnel_size], "Dead", ] <- v_p_SD_tunnels[n_tunnel_size]

# From Dead
a_P["Dead", "Dead", ] <- 1

#******************************************************************************
#### 04.3 Check if transition array and probabilities are valid ####
#******************************************************************************

# Check if transition probability array is valid (i.e., elements cannot < 0 or > 1)
check_transition_probability(a_P,
                             verbose = TRUE) # Prints a message if there's an error

# Check if transition probability array sum to 1 (i.e., each row should sum to 1)
check_sum_of_transition_array(a_P = a_P,
                              n_states = n_states_tunnels,
                              n_cycles = n_t,
                              verbose  = TRUE) # Prints a message if there's an error

#******************************************************************************
#### 05 Run Markov model ####
#******************************************************************************

# Loop through the number of cycles
for (t in 1:n_t) {                     
  # Estimate the Markov trace for cycle t + 1 using the t-th matrix from the 
  # probability array 
  m_M[t + 1, ] <- m_M[t, ] %*% a_P[, , t]  
}

head(m_M, n = 30)

# Create aggregated trace
m_M_tunnels <- cbind(Healthy = m_M[, "Healthy"], 
                     Sick = rowSums(m_M[, 2:(n_tunnel_size + 1)]), 
                     Dead = m_M[, "Dead"])

# Show the first rows of the aggregated Markov trace
head(m_M_tunnels)

# rowSums(m_M_tunnels)

#******************************************************************************
#### 06 Compute and Plot Epidemiological Outcomes ####
#******************************************************************************
#******************************************************************************
#### 06.1 Cohort trace ####
#******************************************************************************

# Create a plot of the data
matplot(m_M_tunnels, 
        type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace", 
        lwd  = 2)              
# Add a legend to the graph
legend("right", 
       v_n, 
       col = c("black", "red", "green"), 
       lty = 1:3, 
       bty = "n")  
# Plot a vertical line that helps identifying at which cycle the prevalence of
# sick is highest.
abline(v = which.max(m_M_tunnels[, "Sick"]) - 1,
       col = "gray")

#******************************************************************************
#### 06.2 Overall Survival (OS) ####
#******************************************************************************

# Calculate the overall survival (OS) probability
v_os <- 1 - m_M_tunnels[, "Dead"]  

# Alternative way of calculating the OS probability  
v_os <- rowSums(m_M_tunnels[, 1:2])   

# Create a simple plot showing the OS
plot(v_os, 
     type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")             
# Add grid 
grid(nx  = n_t, 
     ny  = 10, 
     col = "lightgray", 
     lty = "dotted",
     lwd = par("lwd"),
     equilogs = TRUE) 

#******************************************************************************
#### 06.2.1 Life Expectancy (LE) ####
#******************************************************************************

# Summing probability of OS over time  (i.e. life expectancy)
v_le <- sum(v_os)           

#******************************************************************************
#### 06.3 Disease prevalence ####
#******************************************************************************

# Calculate the disease prevalence
v_prev <- m_M_tunnels[, "Sick"]/v_os

# Plot the disease prevalence
plot(v_prev,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

#******************************************************************************
#### 07 Compute Cost-Effectiveness Outcomes ####
#******************************************************************************
#******************************************************************************
#### 07.1 Mean Costs and QALYs ####
#******************************************************************************

# Per cycle
# Calculate expected costs by multiplying m_M with the cost vector for the different
# Health states
v_tc <- m_M_tunnels %*% v_costs

# Calculate expected QALYs by multiplying m_M with the utilities for the different
# Health states
v_tu <- m_M_tunnels %*% v_utilities

#******************************************************************************
#### 07.2 Discounted Mean Costs and QALYs ####
#******************************************************************************

# Discount costs  by multiplying the cost vector with discount weights (v_dw) 
v_tc_d <-  t(v_tc) %*% v_dwc

# Discount QALYS  by multiplying the QALYs vector with discount weights (v_dw)
v_te_d <-  t(v_tu) %*% v_dwe

#******************************************************************************
#### 07.3 Results ####
#******************************************************************************

results <- data.frame( "Total Discounted Cost" = v_tc_d, 
                       "Life Expectancy" = v_le, 
                       "Total Discounted QALYs" = v_te_d, 
                       check.names = F)
results
