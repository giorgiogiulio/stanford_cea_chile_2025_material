# *****************************************************************************
#
#
# Purpose: Define the 3-State Cancer Relative Survival (CRS) Markov Model
#          as a function for use in calibration exercises
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
# CRS Markov Model Function ----------------------------------------------------
# ******************************************************************************
run_crs_markov <- function(v_params) {
  with(as.list(v_params), {
    ## Markov model parameters
    n_t  <- 60                        # time horizon, number of cycles
    v_n  <- c("NED", "Mets", "Death") # the 3 states of the model
    n_s <- length(v_n)                # number of health states 
    
    # Transition probabilities 
    # p_Mets    = 0.10         	# probability to become sicker when sick
    # p_DieMets = 0.05        	  # hazard ratio of death in sick vs healthy
    
    ####### INITIALIZATION ##########################################
    # create the cohort trace
    m_M <- matrix(NA, nrow = n_t + 1 , 
                  ncol = n_s,
                  dimnames = list(0:n_t, v_n)) # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    
    m_M[1, ] <- c(1, 0, 0)                     # initialize Markov trace
    
    # create transition probability matrix for NO treatment
    m_P <- matrix(0,
                  nrow = n_s, 
                  ncol = n_s,
                  dimnames = list(v_n, v_n))
    # fill in the transition probability array
    ### From NED
    m_P["NED", "NED"]   <- 1 - (p_Mets)
    m_P["NED", "Mets"]  <- p_Mets
    m_P["NED", "Death"] <- 0            # Not allowed to die from cancer in NED state
    ### From Mets
    m_P["Mets", "NED"]   <- 0
    m_P["Mets", "Mets"]  <- 1 - (p_DieMets)
    m_P["Mets", "Death"] <- p_DieMets
    ### From Death
    m_P["Death", "Death"] <- 1
    
    # check rows add up to 1
    if (!isTRUE(all.equal(as.numeric(rowSums(m_P)), as.numeric(rep(1, n_s))))) {
      stop("This is not a valid transition Matrix")
    }
    
    ############# PROCESS ###########################################
    
    for (t in 1:n_t) {                   # throughout the number of cycles
      m_M[t + 1, ] <- m_M[t, ] %*% m_P  # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    
    ####### EPIDEMIOLOGICAL OUTPUT  ###########################################
    #### Overall Survival (OS) ####
    v_os   <- 1 - m_M[, "Death"]  # calculate the overall survival (OS) probability
    #### Proportion of NED and Mets over time ####
    v_ned  <- m_M[, "NED"]
    v_mets <- m_M[, "Mets"]
    #### Ratio of NED-to-Mets over time ####
    r_ned_to_mets <- v_ned/v_mets
    ####### RETURN OUTPUT  ###########################################
    out <- list(Surv = v_os[-c(1:2)]#, r_ned_to_mets = r_ned_to_mets[-c(1:2)]
                )
    
    return(out)
  }
  )
}