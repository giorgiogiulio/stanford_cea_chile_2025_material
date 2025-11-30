# *****************************************************************************
#
# Script: intro_to_r_health_economics.R
#
# Purpose: Introduce students to R and RStudio through hands-on demonstration
#          of basic operations, package installation, data exploration, and
#          cost-effectiveness analysis with health economics examples.
#
# Authors: 
# This work is developed for the Decision Sciences Methods course at:
# Stanford University and Universidad de Chile
#
# - Fernando Alarid-Escudero, PhD
# - Jeremy D. Goldhaber-Fiebert, PhD
# - Jorge Roa, MSc
#
# Course: Laboratorio en R: Proyecto y Costo-efectividad incremental
# Module: Introduction to R and RStudio
#
# *****************************************************************************
#
# Notes:
#
# This script accompanies the "Introduction to RStudio" presentation and is
# designed to be run line-by-line by students to learn R fundamentals using
# health economics and cost-effectiveness examples.
#
# The script demonstrates:
# - Basic R operations and syntax
# - Installing and loading packages (tidyverse, dampack)
# - Creating and exploring health outcome datasets
# - Data manipulation with dplyr
# - Data visualization with ggplot2
# - Cost-effectiveness analysis with dampack
# - Saving outputs following project organization best practices
#
# Students should follow along by:
# 1. Creating a new R Project
# 2. Setting up folder structure (data/, figures/, analysis/)
# 3. Saving this script in the analysis/ folder
# 4. Running code line-by-line using Cmd+Return (Mac) or Ctrl+Enter (PC)
#
# *****************************************************************************

# ******************************************************************************
# 01 Basic R Operations --------------------------------------------------------
# ******************************************************************************

### 01.01 Simple Calculations --------------------------------------------------
# R can be used as a calculator
# Run each line to see the result in the Console

2 + 2
10 * 5
100 / 4

# ******************************************************************************
# 02 Installing and Loading Packages -------------------------------------------
# ******************************************************************************

### 02.01 Installing Packages --------------------------------------------------
# In this section we install all the R packages needed for:
# - Sampling and calibration (lhs, IMIS, matrixStats)
# - Visualization (plotrix, psych, scatterplot3d, ggplot2, GGally)
# - Data manipulation (dplyr)
# - Development tools (devtools, used to install IMIS from the archive)
#
# You only need to install packages **once** on your computer.
# After that, you can comment these lines out.

## Install all required CRAN packages in one shot
# install.packages(c(
#   "dampack",        # cost-effectiveness analysis package
#   "darthtools",     # DARTH tools package
#   "lhs",           # Latin Hypercube Sampling
#   "devtools",      # to install IMIS from archive
#   "matrixStats",   # summary statistics
#   "plotrix",       # for plotCI function
#   "psych",         # for pairs.panels function
#   "scatterplot3d", # 3D visualization
#   "ggplot2",       # general plotting
#   "GGally",        # ggplot-type correlation plots
#   "dplyr"          # data manipulation
# ))


# Use pacman to install (if needed) and load all required packages
pacman::p_load(
  dampack,        # cost-effectiveness analysis package
  darthtools,     # DARTH tools package
  lhs,            # Latin Hypercube Sampling
  devtools,       # development tools; used to install IMIS from archive
  matrixStats,    # fast summary statistics on matrices
  plotrix,        # for plotCI function and extra plotting tools
  psych,          # for pairs.panels and descriptive analyses
  scatterplot3d,  # 3D visualization (scatter plots)
  ggplot2,        # general plotting (grammar of graphics)
  GGally,         # ggplot2 extensions (correlation plots, pairs, etc.)
  dplyr           # data manipulation (filter, mutate, summarise, etc.)
)

# Install IMIS from CRAN archive (only if not already installed)
devtools::install_version( "IMIS", version = "0.1", repos = "http://cran.us.r-project.org" )

### 02.02 Loading Packages -----------------------------------------------------
# Load the packages (you need to do this **every time** you start a new R session)

library(lhs)          # package for Latin Hypercube Sampling
library(IMIS)         # package for Incremental Mixture Importance Sampling
library(matrixStats)  # package used for summary statistics

# visualization
library(plotrix)       # for plotCI function
library(psych)         # for pairs.panels function
library(scatterplot3d) # higher-dimension visualization (3D scatterplots)
library(ggplot2)       # general plotting
library(GGally)        # ggplot-type correlation plots

# data manipulation
library(dplyr)         # for data manipulation (filter, mutate, etc.)

### 02.02 Loading Packages -----------------------------------------------------
# Load the packages (do this every session)

library(tidyverse)
library(dampack)
library(tidyverse)
# You should see messages about the packages loaded


# ******************************************************************************
# 03 Working with Data ---------------------------------------------------------
# ******************************************************************************

### 03.01 Chilean life tables --------------------------------------------

#Load the data

df_lifetable_chile <- read_excel("data/lifetables_chile.xlsx")

### 03.02 Exploring Data Structure ---------------------------------------------
# Always examine your data before analysis

# View entire dataset in a spreadsheet-like window
View(df_lifetable_chile)

# Check the structure of the data
str(df_lifetable_chile)

# Get summary statistics for all variables
summary(df_lifetable_chile)

# View first few rows
head(df_lifetable_chile)


# ******************************************************************************
# 04 Data Exploration with dplyr -----------------------------------------------
# ******************************************************************************

### 04.01 Filtering and Summarizing --------------------------------------------
# Filter students who studied more than 7 hours

# Load required libraries


# Create the data frame (assuming df_lifetable_chile is already loaded)
# If not, you would load it from your data source

# Clean and process the data
df_lifetable_chile_f <- df_lifetable_chile %>%
  filter(`Edad del fallecido` != "Total" & 
           `Edad del fallecido` != "No especificado") %>%
  mutate(
    # Extract numeric age - everything less than 1 year becomes 0
    age_year = case_when(
      str_detect(`Edad del fallecido`, "hora|día|mes") ~ 0,
      str_detect(`Edad del fallecido`, "año") ~ as.numeric(str_extract(`Edad del fallecido`, "\\d+")),
      TRUE ~ NA_real_
    )
  )

# Summary by year
df_lifetable_chile_f_sum <- df_lifetable_chile_f %>%
  group_by(age_year) %>%
  summarise(
    total_deaths = sum(Casos, na.rm = TRUE),
    percentage = (sum(Casos, na.rm = TRUE) / sum(df_lifetable_chile_f$Casos, na.rm = TRUE)) * 100,
    .groups = "drop"
  ) %>%
  arrange(age_year) %>%
  mutate(
    # Create age variable capped at 100
    age = ifelse(age_year > 100, 100, age_year),
    # Create 5-year age groups
    age_group_5yr = case_when(
      age >= 0 & age < 5 ~ "0-4",
      age >= 5 & age < 10 ~ "5-9",
      age >= 10 & age < 15 ~ "10-14",
      age >= 15 & age < 20 ~ "15-19",
      age >= 20 & age < 25 ~ "20-24",
      age >= 25 & age < 30 ~ "25-29",
      age >= 30 & age < 35 ~ "30-34",
      age >= 35 & age < 40 ~ "35-39",
      age >= 40 & age < 45 ~ "40-44",
      age >= 45 & age < 50 ~ "45-49",
      age >= 50 & age < 55 ~ "50-54",
      age >= 55 & age < 60 ~ "55-59",
      age >= 60 & age < 65 ~ "60-64",
      age >= 65 & age < 70 ~ "65-69",
      age >= 70 & age < 75 ~ "70-74",
      age >= 75 & age < 80 ~ "75-79",
      age >= 80 & age < 85 ~ "80-84",
      age >= 85 & age < 90 ~ "85-89",
      age >= 90 & age < 95 ~ "90-94",
      age >= 95 & age <= 100 ~ "95-100"
    )
  )

# Summary by 5-year age groups
df_lifetable_chile_f_sum_5y <- df_lifetable_chile_f_sum%>%
  group_by(age_group_5yr) %>%
  summarise(
    total_deaths = sum(total_deaths),
    percentage = sum(percentage),
    .groups = "drop"
  ) %>%
  mutate(
    # Factorize to maintain logical order
    age_group_5yr = factor(age_group_5yr, levels = c(
      "0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
      "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79",
      "80-84", "85-89", "90-94", "95-100"
    ))
  )


# ******************************************************************************
# 05 Data Visualization with ggplot2 -------------------------------------------
# ******************************************************************************

### 05.01 Scatter Plot ---------------------------------------------------------
# Visualize relationship between hours studied and exam scores

ggplot(df_lifetable_chile_f_sum %>% filter(!is.na(age)), 
             aes(x = age, y = total_deaths)) +
  geom_line(color = "#2E86AB", size = 1) +
  geom_point(color = "#2E86AB", size = 2, alpha = 0.6) +
  labs(
    title = "Life Table: Distribution of Deaths by Age in Chile",
    subtitle = "Ages 0 to 100 (100+ grouped together)",
    x = "Age (years)",
    y = "Number of Deaths",
    caption = paste0("Total deaths: ", format(sum(df_lifetable_chile_f_sum$Casos, na.rm = TRUE), big.mark = ","),
                     "\nNote: Age 0 includes all deaths under 1 year (hours, days, months)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9, color = "gray40")
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10))

# Bar chart by 5-year age groups
ggplot(df_lifetable_chile_f_sum_5y, 
       aes(x = age_group_5yr, y = total_deaths, fill = total_deaths)) +
  geom_col() +
  geom_text(aes(label = paste0(format(total_deaths, big.mark = ","), "\n", 
                               round(percentage, 1), "%")), 
            vjust = -0.3, size = 3) +
  scale_fill_gradient(low = "#A8DADC", high = "#1D3557") +
  labs(
    title = "Deaths by 5-Year Age Groups in Chile",
    x = "Age Group",
    y = "Number of Deaths"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


# ******************************************************************************
# 07 Cost-Effectiveness Analysis with dampack ----------------------------------
# ******************************************************************************

library(dampack)

### 07.01 Example 1: HIV Screening Strategies ----------------------------------
# From: Paltiel et al. 2006. Annals of Internal Medicine.
# Evaluating HIV screening frequencies in high-risk population

## Define strategies
v_hiv_strat_names <- c("Status Quo", "One time", "5 year", "3 year", "1 year")

## Average cost per person for each strategy
v_hiv_costs <- c(26000, 27000, 28020, 28440, 29440)

## Quality-adjusted life-years (QALYs) for each strategy
# Original data in quality-adjusted life-months, converted to years by dividing by 12
v_hiv_qalys <- c(277.25, 277.57, 277.78, 277.83, 277.76) / 12

## Calculate ICERs using dampack
icer_hiv <- dampack::calculate_icers(cost = v_hiv_costs,
                            effect = v_hiv_qalys,
                            strategies = v_hiv_strat_names)

# View results
print(icer_hiv)

#    Strategy  Cost   Effect Inc_Cost  Inc_Effect      ICER Status
# 1 Status Quo 26000 23.10417       NA          NA        NA     ND
# 2   One time 27000 23.13083     1000 0.026666667  37500.00     ND
# 3     5 year 28020 23.14833     1020 0.017500000  58285.71     ND
# 4     3 year 28440 23.15250      420 0.004166667 100800.00     ND
# 5     1 year 29440 23.14667       NA          NA        NA      D

# ------------------------------------------------------------------------------
# How incremental values are calculated
# ------------------------------------------------------------------------------
# The strategies are first ordered by:
#  1) Increasing cost
#  2) Removal of dominated / extendedly dominated strategies
#
# Then, for each non-dominated strategy i (after the first on the frontier):
#   Inc_Cost_i   = Cost_i   - Cost_(i-1)
#   Inc_Effect_i = Effect_i - Effect_(i-1)
#
# where (i-1) is the previous NON-DOMINATED strategy on the cost-effectiveness
# frontier. The ICER is then:
#   ICER_i = Inc_Cost_i / Inc_Effect_i
#
# For the reference (cheapest) strategy and for dominated strategies,
# Inc_Cost, Inc_Effect, and ICER are set to NA.

# ------------------------------------------------------------------------------
# Interpretation of ICER results
# ------------------------------------------------------------------------------
#
# Columns:
# - "Cost"   = expected total cost per person for each strategy.
# - "Effect" = expected health outcome (e.g., QALYs) per person.
# - "Inc_Cost"   = additional cost vs. the previous non-dominated strategy.
# - "Inc_Effect" = additional QALYs vs. the previous non-dominated strategy.
# - "ICER"   = incremental cost-effectiveness ratio = Inc_Cost / Inc_Effect.
# - "Status" = ND = non-dominated; D = dominated.
#
# Row 1: Status Quo
# - Cheapest option (Cost = 26,000) and baseline Effect = 23.10.
# - No incremental values (Inc_Cost, Inc_Effect, ICER = NA) because it is
#   the reference strategy on the frontier.
#
# Row 2: One time
# - Compared with Status Quo:
#   Inc_Cost   = 27,000 - 26,000 = 1,000
#   Inc_Effect = 23.13083 - 23.10417 ≈ 0.02667
# - ICER = 1,000 / 0.02667 ≈ 37,500 per QALY.
# - Marked ND (non-dominated), so it remains on the cost-effectiveness frontier.
#
# Row 3: 5 year
# - Compared with One time:
#   Inc_Cost   = 28,020 - 27,000 = 1,020
#   Inc_Effect = 23.14833 - 23.13083 = 0.01750
# - ICER ≈ 1,020 / 0.01750 ≈ 58,286 per QALY.
# - Still ND, so also on the frontier.
#
# Row 4: 3 year
# - Compared with 5 year:
#   Inc_Cost   = 28,440 - 28,020 = 420
#   Inc_Effect = 23.15250 - 23.14833 ≈ 0.00417
# - ICER ≈ 420 / 0.00417 ≈ 100,800 per QALY.
# - ND, but with a much higher ICER (more expensive per extra QALY).
#
# Row 5: 1 year
# - Costs more (29,440) but has LOWER Effect (23.14667) than the 3 year
#   strategy (23.15250).
# - It is strictly dominated (more costly and less effective),
#   so Status = "D" and Inc_Cost / Inc_Effect / ICER are NA.


### 07.02 Visualize HIV CEA Results --------------------------------------------

# Basic cost-effectiveness plane
plot(icer_hiv)

# Basic cost-effectiveness plane and flip the axes for better visualization
plot(icer_hiv) +
  coord_flip()

# Plot with all strategies labeled. it's not by default. 
plot(icer_hiv, label = "all")

# Customized plot with title
plot(icer_hiv, label = "all") +
  theme_classic() +
  ggtitle("Cost-effectiveness of HIV screening strategies")


### 07.03 Example 2: C. difficile Treatment Strategies ------------------------
# From: Rajasingham et al. 2020. Clinical Infectious Diseases.
# Using probabilistic sensitivity analysis data

## Load C. diff PSA data (included in dampack package)
data("psa_cdiff")


## Calculate mean cost and effectiveness for each strategy
df_cdiff_ce <- dampack:::summary.psa(psa_cdiff,)

head(df_cdiff_ce)

## Calculate ICERs
icer_cdiff <- calculate_icers(cost = df_cdiff_ce$meanCost,
                              effect = df_cdiff_ce$meanEffect,
                              strategies = df_cdiff_ce$Strategy)



# View all results
print(icer_cdiff)

plot(icer_cdiff, label = "all") +
  coord_flip()

# View only non-dominated strategies
icer_cdiff %>%
  filter(Status == "ND")


### 07.04 Visualize C. diff CEA Results ----------------------------------------

# Basic plot
plot(icer_cdiff)

# Plot with all strategies labeled
plot(icer_cdiff, label = "all")


# Plot only efficient frontier (non-dominated strategies)
plot(icer_cdiff, plot_frontier_only = TRUE)

# Customized plot with axis labels
plot(icer_cdiff, 
     currency = "USD", 
     effect_units = "quality-adjusted life-years")

#   Strategy     Cost   Effect Inc_Cost Inc_Effect      ICER Status
# 1        s3 57336.01 12.93996       NA         NA        NA     ND
# 2       s27 57541.25 13.01406 205.2466 0.07410015  2769.855     ND
# 3       s33 57642.26 13.03891 101.0061 0.02484756  4065.031     ND
# 4       s31 57934.07 13.09663 291.8156 0.05771416  5056.222     ND
# 5       s43 58072.11 13.11286 138.0394 0.01623188  8504.216     ND
# 6       s44 58665.78 13.12833 593.6686 0.01547517 38362.652     ND
# 7       s39 57814.65 13.04628       NA         NA        NA     ED
# 8        s4 57887.48 12.99707       NA         NA        NA      D
# 9       s13 58018.63 13.06504       NA         NA        NA      D
#10       s37 58081.79 13.10297       NA         NA        NA      D
#11      s20 58634.20 13.11006       NA         NA        NA      D

# ------------------------------------------------------------------------------
# Interpretation of PSA ICER table
# ------------------------------------------------------------------------------
# This table summarizes results from the probabilistic sensitivity analysis (PSA).
# Each row represents the *expected* (mean) cost and effect across PSA iterations.
#
# The strategies are ordered by increasing expected cost, and then classified as:
# - ND: Non-dominated (on the cost-effectiveness frontier)
# - ED: Extendedly dominated
# - D : Dominated

# Non-dominated strategies (ND)
# Strategies s3, s27, s33, s31, s43, and s44 form the cost-effectiveness frontier.
# Moving along the frontier:
# - Costs and effects both increase
# - ICERs increase monotonically, as expected under efficiency

# Example interpretation:
# - s27 is the first improvement over s3 and has a very low ICER (~2,770),
#   meaning large health gains for a small increase in cost.
# - s33 and s31 provide additional gains at higher ICERs.
# - s44 is the most effective strategy but also much more expensive,
#   with an ICER of ~38,000 per unit of effect.

# Extended dominance (ED)
# Strategy s39 is extendedly dominated:
# - It is less efficient than a linear combination of two ND strategies
# - Even though it improves health, it is excluded from the optimal frontier

#  Dominated (D)
# Strategies s4, s13, s37, and s20 are strictly dominated:
# - They cost more and produce fewer effects than at least one other strategy
# - They are never optimal at any willingness-to-pay threshold


# Along the cost-effectiveness frontier, ICERs must be monotonically increasing.
# This means that each more effective strategy should have a higher (or equal)
# ICER than the previous one.
#
# If an ICER decreases after increasing (i.e., goes up and then down),
# the intermediate strategy must be removed from the frontier.
# Such a strategy is said to be *extendedly dominated* and is excluded
# from cost-effectiveness consideration.

# ******************************************************************************
# 09 Key Takeaways -------------------------------------------------------------
# ******************************************************************************

### 09.01 R Programming Fundamentals -------------------------------------------
# 1. Use <- to assign values to objects
# 2. Objects appear in the Environment panel
# 3. Functions are called with parentheses: function(arguments)
# 4. Comments start with # and are ignored by R

### 09.02 Package Management ---------------------------------------------------
# 5. Install packages ONCE with install.packages("package_name")
# 6. Load packages EVERY SESSION with library(package_name)
# 7. tidyverse includes many useful packages for data analysis
# 8. dampack provides specialized functions for cost-effectiveness analysis

### 09.03 Data Manipulation ----------------------------------------------------
# 9. Use the pipe operator %>% to chain operations
# 10. dplyr provides intuitive functions: filter(), mutate(), group_by(), summarise()
# 11. Always explore your data before analysis: str(), summary(), View()

### 09.04 Visualization --------------------------------------------------------
# 12. ggplot2 creates beautiful, publication-ready visualizations
# 13. Basic template: ggplot(data, aes(x, y)) + geom_*() + labs() + theme_*()
# 14. Save plots with ggsave()

### 09.05 Cost-Effectiveness Analysis ------------------------------------------
# 15. Use calculate_icers() from dampack to conduct CEA
# 16. ICERs compare incremental costs to incremental effects
# 17. Dominated strategies are identified automatically (D = strong, ED = weak)
# 18. plot() method creates cost-effectiveness planes
# 19. Always consider all feasible strategies in your analysis

### 09.06 Project Organization -------------------------------------------------
# 20. Follow a consistent folder structure (data/, figures/, analysis/)
# 21. Save processed data with write.csv()
# 22. Document your work with comments


# ******************************************************************************
# 10 Next Steps ----------------------------------------------------------------
# ******************************************************************************

# CONGRATULATIONS! You've completed the introduction to R!
#
# To practice further:
# 1. Modify the student dataset (add more students, variables)
# 2. Try the HIV example with different costs or QALYs
# 3. Explore the full C. diff PSA data with psa_cdiff
# 4. Load real datasets from CSV files
# 5. Create your own cost-effectiveness analysis
#
# Resources for continued learning:
# - R for Data Science: https://r4ds.had.co.nz/
# - DARTH group materials: http://darthworkgroup.com/
# - dampack vignettes: vignette("basic_cea", package = "dampack")
# - RStudio Cheatsheets: https://www.rstudio.com/resources/cheatsheets/
#
# *****************************************************************************
# END OF SCRIPT
# *****************************************************************************