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
### MARKOV MODEL - THREE STATE ####
#******************************************************************************

# Clear memory (removes all the variables from the workspace)
rm(list = ls()) 

#******************************************************************************
#### 01 Load packages ####
#******************************************************************************

# Load packages from CRAN
library("diagram") 

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
v_n  <- c("Healthy", "Sick", "Dead")    # State names
n_states  <- length(v_n)                # Number of states
n_age_init <- 25  # age at baseline
n_age_max  <- 100 # maximum age of follow up
n_t        <- n_age_max - n_age_init # Number of cycles

p_HD <- 0.05                            # Probability to die when healthy
p_HS <- 0.20                            # Probability to become sick when healthy
p_SD <- 0.15                            # Probability to die when sick

# Costs and utilities  
c_H         <- 400                      # Cost of remaining one cycle healthy
c_S         <- 1000                     # Cost of remaining one cycle sick
c_D         <- 0                        # Cost of remaining one cycle dead
u_H         <- 1                        # Utility when healthy 
u_S         <- 0.5                      # Utility when sick
u_D         <- 0                        # Utility when dead
d_e         <- d_c <- 0.03              # Equal discount of costs and QALYs by 3%
v_costs     <- c(c_H, c_S, c_D)         # All costs
v_utilities <- c(u_H, u_S, u_D)         # All utilities

# Calculate discount weights for costs for each cycle based on discount rate d_c
v_dwc <- 1 / (1 + d_e) ^ (0:n_t) 

# Calculate discount weights for effectiveness for each cycle based on discount rate d_e
v_dwe <- 1 / (1 + d_c) ^ (0:n_t) 

## Draw the state-transition cohort model
m_P_diag <- matrix(0, 
                   nrow     = n_states, 
                   ncol     = n_states, 
                   dimnames = list(v_n, v_n))
m_P_diag["Healthy", "Sick" ]     = "" 
m_P_diag["Healthy", "Dead" ]     = ""
m_P_diag["Healthy", "Healthy" ]  = ""
m_P_diag["Sick"   , "Dead" ]     = ""
m_P_diag["Sick"   , "Sick" ]     = ""
m_P_diag["Dead"   , "Dead" ]     = ""
layout.fig <- c(2, 1)
plotmat(A        = t(m_P_diag), 
        pos      = t(layout.fig), 
        self.cex = 0.5, 
        curve    = 0, 
        arr.pos  = 0.8,  
        latex    = T,
        arr.type = "curved", 
        relsize  = 0.85, 
        box.prop = 0.8, 
        cex      = 0.8, 
        box.cex  = 0.7, 
        lwd      = 1)

#******************************************************************************
#### 04 Define and initialize matrices and vectors ####
#******************************************************************************
#******************************************************************************
#### 04.1 Cohort trace ####
#******************************************************************************

# Create the cohort trace
m_M <- matrix(NA, 
              # Create Markov trace (n_t + 1 because R doesn't understand Cycle 0)
              nrow = n_t + 1,
              ncol = n_states, 
              dimnames = list(0:n_t, v_n))

# Initialize first cycle of Markov trace
m_M[1, ] <- c(1, 0, 0)

#******************************************************************************
#### 04.2 Transition probability matrix ####
#******************************************************************************

# Create the transition probability matrix
m_P  <- matrix(0,
               nrow     = n_states, 
               ncol     = n_states,
               dimnames = list(v_n, v_n)) # Name the columns and rows of the transition 

# Probability matrix
m_P

#Fill in the transition probability matrix:
# From Healthy
m_P["Healthy", "Healthy"] <- 1 - p_HD - p_HS
m_P["Healthy", "Sick"]    <- p_HS
m_P["Healthy", "Dead"]    <- p_HD

# From Sick
m_P["Sick", "Sick"] <- 1 - p_SD
m_P["Sick", "Dead"] <- p_SD

# From Dead
m_P["Dead", "Dead"] <- 1

# Check rows add up to 1
rowSums(m_P)

#******************************************************************************
#### 05 Run Markov model ####
#******************************************************************************

# Loop through the number of cycles
for (t in 1:n_t) {
  # Estimate the state vector for the next cycle (t + 1)
  m_M[t + 1, ] <- m_M[t, ] %*% m_P  
}

#******************************************************************************
#### 06 Compute and Plot Epidemiological Outcomes ####
#******************************************************************************
#******************************************************************************
#### 06.1 Cohort trace ####
#******************************************************************************

# Create a plot of the data
matplot(x    = c(0:n_t),
        y    = m_M, 
        type = 'l',
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Cohort Trace", 
        lwd  = 3)
# Add a legend to the graph
legend("right",                                 
       legend = v_n, 
       col    = c("black", "red", "green"), 
       lty    = 1:3, 
       bty    = "n")                            
# Plot a vertical line that helps identifying at which cycle the prevalence of 
# sick is highest.  
abline(v = which.max(m_M[, "Sick"]) - 1, 
       col = "gray") 

#******************************************************************************
#### 06.2 Overall Survival (OS) ####
#******************************************************************************

# Calculate the overall survival (OS) probability
v_os <- 1 - m_M[, "Dead"]

# Alternative way of calculating the OS probability
v_os <- rowSums(m_M[, 1:2])

# Create a simple plot showing the OS
plot(x    = c(0:n_t),
     y    = v_os,
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

# Summing probability of OS over time (i.e. life expectancy)
v_le <- sum(v_os)           

#******************************************************************************
#### 06.3 Disease prevalence ####
#******************************************************************************

# Calculate the disease prevalence
v_prev <- m_M[, "Sick"]/v_os

# Plot the disease prevalence
plot(x    = c(0:n_t), 
     y    = v_prev,
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
v_tc <- m_M %*% v_costs

# Calculate expected QALYs by multiplying m_M with the utilities for the different 
# Health states   
v_tu <- m_M %*% v_utilities  

#******************************************************************************
#### 07.2 Discounted Mean Costs and QALYs ####
#******************************************************************************

# Discount costs by multiplying the cost vector with discount weights (v_dw)
n_tc_d <-  t(v_tc) %*% v_dwc
#sum(v_tc * v_dwc)

# Discount QALYS by multiplying the QALYs vector with discount weights (v_dw)
n_te_d <-  t(v_tu) %*% v_dwe
#sum(v_tu * v_dwe)

#******************************************************************************
#### 07.3 Results ####
#******************************************************************************

results <- data.frame( "Total Discounted Cost" = n_tc_d, 
                       "Life Expectancy" = v_le, 
                       "Total Discounted QALYs" = n_te_d, 
                       check.names = F)
results
