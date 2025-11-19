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

### 01.02 Creating Objects (Variables) -----------------------------------------
# Use the assignment operator <- to store values
# These objects will appear in your Environment panel

my_number <- 42
my_name <- "Student"
my_vector <- c(1, 2, 3, 4, 5)

# View objects in Environment panel
print(my_number)
print(my_name)
print(my_vector)

### 01.03 Basic Operations with Vectors ----------------------------------------
# Vectors are the basic data structure in R

sum(my_vector)
mean(my_vector)
max(my_vector)


# ******************************************************************************
# 02 Installing and Loading Packages -------------------------------------------
# ******************************************************************************

### 02.01 Installing Packages --------------------------------------------------
# Install tidyverse and dampack (only need to do this once)
# Uncomment the lines below if you haven't installed these packages yet

# install.packages("tidyverse")
# install.packages("dampack")

### 02.02 Loading Packages -----------------------------------------------------
# Load the packages (do this every session)

library(tidyverse)
library(dampack)

# You should see messages about the packages loaded


# ******************************************************************************
# 03 Working with Data ---------------------------------------------------------
# ******************************************************************************

### 03.01 Creating a Sample Dataset --------------------------------------------
# We'll create a dataset of student study habits and exam performance

student_data <- data.frame(
  student_id = 1:20,
  age = c(22, 25, 23, 24, 22, 26, 23, 25, 24, 22,
          23, 24, 25, 22, 26, 23, 24, 25, 22, 23),
  hours_studied = c(5, 8, 6, 7, 4, 9, 6, 8, 7, 5,
                    6, 7, 8, 5, 10, 6, 7, 8, 4, 6),
  exam_score = c(65, 85, 70, 78, 60, 92, 72, 88, 80, 68,
                 73, 79, 86, 67, 95, 74, 81, 87, 62, 75)
)

### 03.02 Exploring Data Structure ---------------------------------------------
# Always examine your data before analysis

# View entire dataset in a spreadsheet-like window
View(student_data)

# Check the structure of the data
str(student_data)

# Get summary statistics for all variables
summary(student_data)

# View first few rows
head(student_data)


# ******************************************************************************
# 04 Data Exploration with dplyr -----------------------------------------------
# ******************************************************************************

### 04.01 Filtering and Summarizing --------------------------------------------
# Filter students who studied more than 7 hours

high_study <- student_data %>%
  filter(hours_studied > 7)

print(high_study)

# Calculate average exam score by study hours category
study_summary <- student_data %>%
  mutate(study_category = case_when(
    hours_studied < 6 ~ "Low",
    hours_studied >= 6 & hours_studied < 8 ~ "Medium",
    hours_studied >= 8 ~ "High"
  )) %>%
  group_by(study_category) %>%
  summarise(
    avg_score = mean(exam_score),
    count = n()
  )

print(study_summary)


# ******************************************************************************
# 05 Data Visualization with ggplot2 -------------------------------------------
# ******************************************************************************

### 05.01 Scatter Plot ---------------------------------------------------------
# Visualize relationship between hours studied and exam scores

ggplot(student_data, aes(x = hours_studied, y = exam_score)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    title = "Relationship between Study Hours and Exam Scores",
    x = "Hours Studied",
    y = "Exam Score"
  ) +
  theme_minimal()

### 05.02 Histogram ------------------------------------------------------------
# Distribution of exam scores

ggplot(student_data, aes(x = exam_score)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  labs(
    title = "Distribution of Exam Scores",
    x = "Exam Score",
    y = "Frequency"
  ) +
  theme_minimal()

### 05.03 Box Plot -------------------------------------------------------------
# Exam scores by age

ggplot(student_data, aes(x = factor(age), y = exam_score)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(
    title = "Exam Scores by Age",
    x = "Age",
    y = "Exam Score"
  ) +
  theme_minimal()


# ******************************************************************************
# 06 Statistical Analysis ------------------------------------------------------
# ******************************************************************************

### 06.01 Correlation Analysis -------------------------------------------------
# Examine the relationship between hours studied and exam score

correlation <- cor(student_data$hours_studied, student_data$exam_score)
print(paste("Correlation coefficient:", round(correlation, 3)))

### 06.02 Linear Regression Model ----------------------------------------------
# Build a model to predict exam scores from study hours

model <- lm(exam_score ~ hours_studied, data = student_data)
summary(model)

# Visualize the regression
ggplot(student_data, aes(x = hours_studied, y = exam_score)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") +
  labs(
    title = "Linear Regression: Study Hours Predict Exam Scores",
    subtitle = paste("R² =", round(summary(model)$r.squared, 3)),
    x = "Hours Studied",
    y = "Exam Score"
  ) +
  theme_minimal()


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

dampack:::summary.psa_cdiff

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


# ******************************************************************************
# 08 Saving Your Work ----------------------------------------------------------
# ******************************************************************************

### 08.01 Create Necessary Folders ---------------------------------------------
# Create folders if they don't exist

if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("data")) dir.create("data")

### 08.02 Save Plots -----------------------------------------------------------
# Save student data plots

ggsave("figures/study_score_relationship.png", 
       width = 8, height = 6, dpi = 300)

# Save HIV CEA plot
hiv_plot <- plot(icer_hiv, label = "all") +
  theme_classic() +
  ggtitle("Cost-effectiveness of HIV screening strategies")

ggsave("figures/hiv_cea_plot.png", 
       plot = hiv_plot,
       width = 8, height = 6, dpi = 300)

# Save C. diff CEA plot
cdiff_plot <- plot(icer_cdiff, plot_frontier_only = TRUE)

ggsave("figures/cdiff_cea_plot.png", 
       plot = cdiff_plot,
       width = 8, height = 6, dpi = 300)

print("All plots saved to figures/ folder")

### 08.03 Save Data ------------------------------------------------------------
# Save student data
write.csv(student_data, "data/student_data.csv", row.names = FALSE)
write.csv(study_summary, "data/study_summary.csv", row.names = FALSE)

# Save CEA results
write.csv(icer_hiv, "data/hiv_cea_results.csv", row.names = FALSE)
write.csv(icer_cdiff, "data/cdiff_cea_results.csv", row.names = FALSE)

print("All data saved to data/ folder")


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