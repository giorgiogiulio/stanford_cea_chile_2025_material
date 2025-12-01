#
# 00 Script Metadata and Context -----------------------------------------------
#
# Script: intro_to_r_health_economics.R#
# Purpose: Introduce students to R and RStudio through hands-on demonstration
#          of basic operations, package installation, data exploration, and
#          cost-effectiveness analysis with health economics examples.
#
# This work is developed for the Decision Sciences Methods course at:
# Stanford University and Universidad de Chile
#
# Contributors:
# - Jorge Roa, MSc
# - Jeremy D. Goldhaber-Fiebert, PhD
# - Fernando Alarid-Escudero, PhD
#
# Course: Laboratorio en R: Costo-efectividad incremental
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

#Execute code in R

#Command + Enter (Mac) or Ctrl + Enter (PC) to run each line

# 01 Basic R Operations --------------------------------------------------------

## 01.01 Simple Calculations ---------------------------------------------------
# R can be used as a calculator
# Run each line to see the result in the Console

2 + 2
10 * 5
100 / 4

## 01.02 Creating Objects ------------------------------------------------------

# Assign values to objects using <-

#Command + Shift + M to create the <- operator in Mac
#Ctrl + Shift + M to create the <- operator in PC

a <- 10

b <- 5

# Perform operations using objects

c <- a + b

# View the value of an object by typing its name

c

# 02 Installing and Loading Packages -------------------------------------------

## 02.01 Installing Packages (CRAN) --------------------------------------------
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

install.packages("cli") #For windows support. 
install.packages("rlang") # package management tool
install.packages("devtools") # development tools; used to install IMIS from archive
install.packages("pacman") # development tools; used to install IMIS from archive

install.packages("dampack") #For windows support. 
devtools::install_github("DARTH-git/darthtools")
library(devtools)

# Use pacman to install (if needed) and load all required packages
pacman::p_load(       
  lhs,            # Latin Hypercube Sampling
  matrixStats,    # fast summary statistics on matrices
  plotrix,        # for plotCI function and extra plotting tools
  psych,          # for pairs.panels and descriptive analyses
  scatterplot3d,  # 3D visualization (scatter plots)
  ggplot2,        # general plotting (grammar of graphics)
  GGally,         # ggplot2 extensions (correlation plots, pairs, etc.)
  dplyr,          # data manipulation (filter, mutate, summarise, etc.)
  readxl          # to read and write excel files
)



library(cli)
library(dampack) # cost-effectiveness analysis package
library(darthtools) # DARTH tools package


# Install IMIS from CRAN archive (only if not already installed)
devtools::install_version( "IMIS", version = "0.1", 
                           repos = "http://cran.us.r-project.org" )

## 02.02 Loading Packages (Long List) ------------------------------------------
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

## 02.03 Loading Tidyverse and Dampack -----------------------------------------
# Load the packages (do this every session)

library(tidyverse)
library(dampack)
library(tidyverse)
library(readxl)
# You should see messages about the packages loaded


# 03 Working with Data ---------------------------------------------------------

## 03.01 Chilean life tables ---------------------------------------------------

#Load the data

df_lifetable_chile <- read_excel("data/lifetables_chile.xlsx")

## 03.02 Exploring Data Structure ----------------------------------------------
# Always examine your data before analysis

# View entire dataset in a spreadsheet-like window
View(df_lifetable_chile)

# Check the structure of the data
str(df_lifetable_chile)

# Get summary statistics for all variables
summary(df_lifetable_chile)

# View first 6 rows
head(df_lifetable_chile)


# 04 Data Manipulation with dplyr ----------------------------------------------

## 04.01 Filtering and Summarizing ---------------------------------------------
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
df_lifetable_chile_f_sum_5y <- df_lifetable_chile_f_sum %>%
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



# 05 Data Visualization with ggplot2 -------------------------------------------

## 05.01 Line and Scatter Plot (Yearly Age) ------------------------------------
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
    caption = paste0("Total deaths: ", 
                     format(sum(df_lifetable_chile_f_sum$Casos, na.rm = TRUE), 
                            big.mark = ","),
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

## 05.02 Bar Chart (5-Year Age Groups) -----------------------------------------
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


# 06 Cost-Effectiveness Analysis ----------------------------
# 07 Cost-Effectiveness Analysis with dampack ----------------------------------
library(dampack)

## 07.01 Example 1: Cervical Cancer Screening (Chilean Case) --------------------
# Contexto: En Chile, el cáncer cervicouterino (CaCu) es parte del plan AUGE/GES
# desde 2005. Actualmente, el programa nacional se basa en el Papanicolaou (PAP)
# cada 3 años para mujeres de 25-64 años. Sin embargo, estudios chilenos han
# demostrado que el test de VPH (Virus Papiloma Humano) tiene mayor sensibilidad.
#
# Referencias:
# - Ferreccio et al. (2013): Estudio en Santiago comparando PAP vs VPH
#   * Sensibilidad PAP: 22.1% | VPH: 92.7% para detectar CIN2+
# - MINSAL: ~600 mujeres mueren anualmente de CaCu en Chile
# - Cobertura actual de PAP: <70% de la población objetivo

### 07.01.01 Define Strategies and Data ----------------------------------------
## Definir estrategias de tamizaje
v_cacu_strat_names<- c("Sin tamizaje", 
                       "PAP cada 5 años", 
                       "PAP cada 3 años (actual)", 
                       "VPH cada 5 años",
                       "VPH cada 3 años")

## Costos promedio por mujer durante su vida (en pesos chilenos)
# Basados en costos del sistema público de salud (FONASA)
# Incluyen: costo de examen, seguimiento, colposcopía, tratamiento
v_cacu_costs <- c(50000,      # Sin tamizaje (solo costos de tratamiento tardío)
                  320000,     # PAP cada 5 años
                  450000,     # PAP cada 3 años (estrategia actual)
                  380000,     # VPH cada 5 años (test más caro, menos frecuente)
                  520000)     # VPH cada 3 años

## Años de Vida Ajustados por Calidad (QALYs o QALYs)
# Basados en esperanza de vida en Chile (~80 años para mujeres)
# y calidad de vida asociada a prevenir Cancer Cervicouterino

v_cacu_qalys <- c(22.5,   # Sin tamizaje (mayor mortalidad por CaCu)
                  24.2,   # PAP cada 5 años
                  24.4,   # PAP cada 3 años
                  24.8,   # VPH cada 5 años (mayor sensibilidad)
                  24.9)   # VPH cada 3 años

### 07.01.02 Calculate ICERs and Plot Results ----------------------------------
## Calcular icers usando dampack
icer_cacu <- dampack::calculate_icers(cost = v_cacu_costs,
                                      effect = v_cacu_qalys,
                                      strategies = v_cacu_strat_names)
# Ver resultados
print(icer_cacu)

# Gráfico personalizado con título y tema
plot(icer_cacu, label = "all") +
  theme_minimal() +
  labs(title = "Costo-Efectividad del Tamizaje de Cáncer Cervicouterino en Chile",
       subtitle = "Comparación de estrategias para el plan AUGE/GES",
       x = "QALYs (Años de Vida Ajustados por Calidad)",
       y = "Costo (Pesos Chilenos)",
       caption = "Basado en datos del sistema público de salud (FONASA)")


# Columnas:
# - "Cost"       = costo esperado por mujer durante su vida (pesos chilenos)
# - "Effect"     = QALYs esperados por mujer
# - "Inc_Cost"   = costo incremental vs. estrategia no dominada anterior
# - "Inc_Effect" = QALYs incrementales vs. estrategia no dominada anterior
# - "ICER"       = razón de costo-efectividad incremental (pesos por QALY ganado)
# - "Status"     = ND = no dominada; D = dominada

#_______________________________________________________________________________
#                   Strategy   Cost Effect Inc_Cost Inc_Effect      ICER Status
# 1             Sin tamizaje  50000   22.5       NA         NA        NA     ND
# 2          VPH cada 5 años 380000   24.8   330000        2.3  143478.3     ND
# 3          VPH cada 3 años 520000   24.9   140000        0.1 1400000.0     ND
# 4          PAP cada 5 años 320000   24.2       NA         NA        NA     ED
# 5 PAP cada 3 años (actual) 450000   24.4       NA         NA        NA      D

### 07.01.03 Incremental Calculations Explanation ------------------------------
# Cómo se calculan los valores incrementales
# Las estrategias se ordenan primero por:
#  1) Costo creciente (Increasing cost)
#  2) Eliminación de estrategias dominadas / extendidamente dominadas
#     (Removal of dominated / extendedly dominated strategies)
#
# Luego, para cada estrategia no dominada i
# (exceptuando la primera en la frontera):
#   Inc_Cost_i   = Cost_i   - Cost_(i-1)
#   Inc_Effect_i = Effect_i - Effect_(i-1)
#
# donde (i-1) es la estrategia NO DOMINADA previa en la
# frontera de costo-efectividad (cost-effectiveness frontier).
# El ICER se calcula como:
#   ICER_i = Inc_Cost_i / Inc_Effect_i
#
# Para la estrategia de referencia (la más barata)
# y para las estrategias dominadas,
# Inc_Cost, Inc_Effect e ICER se asignan como NA.

### 07.01.04 Detailed Interpretation of Results --------------------------------
#
# FILA 1: Sin tamizaje (estrategia de referencia)
# - Costo más bajo: $50,000 por mujer
# - Menor efectividad: 22.5 QALYs
# - Status ND: Es la referencia en la frontera de costo-efectividad
# - Interpretación: Aunque es la más barata, implica alta mortalidad por CaCu
#   y menor calidad de vida. No es éticamente aceptable como política pública.
#
# FILA 2: VPH cada 5 años
# - Costo: $380,000 | Efectividad: 24.8 QALYs
# - Inc_Cost: $330,000 adicionales vs. "Sin tamizaje"
# - Inc_Effect: 2.3 QALYs adicionales vs. "Sin tamizaje"
# - icer: $143,478 por QALY ganado
# - Status ND: No dominada, permanece en la frontera
# - Interpretación:
#   * El icer de ~$143,000/QALY está muy por debajo del umbral típico
#     de 1 PIB per cápita (~$10 millones en Chile)
#   * Gana 2.3 QALYs vs. no hacer nada (equivalente a ~2.3 años de vida
#     saludable adicionales por mujer)
#   * El test VPH tiene mayor sensibilidad (92.7% vs 22.1% del PAP)
#   * Permite intervalos de 5 años (vs 3 del PAP) por mayor sensibilidad
#   * Posibilita autotoma → aumenta cobertura en zonas rurales
#
# FILA 3: VPH cada 3 años
# - Costo: $520,000 | Efectividad: 24.9 QALYs
# - Inc_Cost: $140,000 adicionales vs. "VPH cada 5 años"
# - Inc_Effect: 0.1 QALYs adicionales vs. "VPH cada 5 años"
# - icer: $1,400,000 por QALY ganado
# - Status ND: No dominada, pero con icer muy alto
# - Interpretación: Esta estrategia NO es costo-efectiva
#   * Cuesta $1.4 millones adicionales para ganar apenas 0.1 QALYs
#     (equivalente a ~36 días de vida adicionales)
#   * El icer de $1.4M/QALY supera cualquier umbral razonable
#   * La ganancia marginal es mínima comparada con VPH cada 5 años
#   * Desde perspectiva de política pública: NO RECOMENDABLE
#
# FILA 4: PAP cada 5 años
# - Costo: $320,000 | Efectividad: 24.2 QALYs
# - Status ED: Extendidamente dominada
# - Interpretación: ¡La estrategia con PAP está EXTENDIDAMENTE DOMINADA!
#   * Cuesta $320,000 y da 24.2 QALYs
#   * "VPH cada 5 años" cuesta $380,000 y da 24.8 QALYs
#   * Diferencia: $60,000 más para ganar 0.6 QALYs adicionales
#   * icer implícito: $60,000/0.6 = $100,000 por QALY
#   * Este icer es MEJOR que la combinación lineal de otras estrategias
#   * Por lo tanto, PAP cada 5 años es ineficiente (ED)
#
# FILA 5: PAP cada 3 años (actual)
# - Costo: $450,000 | Efectividad: 24.4 QALYs
# - Status D: Estrictamente dominada
# - Interpretación: ¡LA ESTRATEGIA ACTUAL DE CHILE ESTÁ DOMINADA!
#   * "VPH cada 5 años" cuesta MENOS ($380,000 vs $450,000)
#   * "VPH cada 5 años" es MÁS efectivo (24.8 vs 24.4 QALYs)
#   * PAP cada 3 años es más caro Y menos efectivo = DOMINANCIA ESTRICTA
#   * Esto representa una ineficiencia del sistema actual
#   * Chile está gastando $70,000 más por mujer para obtener PEOR resultado
#


### 07.01.05 Policy Implications (Chile) ---------------------------------------
# Interpretación de los resultados para Chile
#
# Implicancias para política pública:
# - Chile podría considerar transición de PAP a VPH como test primario
# - VPH permite intervalos más largos (cada 5 años) manteniendo efectividad
# - Posibilidad de autotoma (mujer toma su propia muestra) aumenta cobertura
# - Especialmente relevante para zonas rurales y población vulnerable


# Notas adicionales sobre el contexto chileno
#
# Desafíos actuales:
# - Cobertura de PAP ha disminuido en última década (<70%)
# - Mortalidad por CaCu en Chile es 4 veces mayor que en países desarrollados
# - Inequidades en acceso: mujeres de zonas rurales y NSE bajo
# - Listas de espera para colposcopía y tratamiento
#
# Oportunidades:
# - Test VPH permite autotoma → mayor alcance poblacional
# - Vacunación VPH para niñas <13 años desde 2015
# - Evidencia local robusta para sustentar cambio de política

# IMPLICANCIAS PARA POLÍTICA PÚBLICA EN CHILE:

# RECOMENDACIÓN PRINCIPAL:
# Chile debería reemplazar "PAP cada 3 años" por "VPH cada 5 años"
#
# Beneficios del cambio:
# 1. AHORRO: $70,000 por mujer ($450k - $380k)
#    - Con ~5 millones de mujeres objetivo → ahorro potencial de $350 mil millones
# 2. EFECTIVIDAD: +0.4 QALYs por mujer (24.8 - 24.4)
#    - Equivalente a ~5 meses adicionales de vida saludable por mujer
# 3. COBERTURA: Autotoma de VPH facilita acceso en zonas rurales
# 4. EQUIDAD: Mayor alcance en población vulnerable y NSE bajo
# 5. LOGÍSTICA: Intervalos de 5 años reducen carga sobre sistema de salud


## 07.02 Example 2: C. difficile Treatment Strategies (PSA) --------------------
# From: Rajasingham et al. 2020. Clinical Infectious Diseases.
# Using probabilistic sensitivity analysis data

### 07.02.01 Load Data and Calculate Mean CE -----------------------------------
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


## 07.03 Visualize C. diff CEA Results -----------------------------------------

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

### 07.03.01 Interpretation of PSA ICER table ----------------------------------
# Interpretation of PSA ICER table
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

# 08 Placeholder Section (PSA Visualization/Analysis) --------------------------
# NOTE: This section number was skipped in the original document (jump from 07 to 09).

# 09 Key Takeaways -------------------------------------------------------------

## 09.01 R Programming Fundamentals --------------------------------------------
# 1. Use <- to assign values to objects
# 2. Objects appear in the Environment panel
# 3. Functions are called with parentheses: function(arguments)
# 4. Comments start with # and are ignored by R

## 09.02 Package Management ----------------------------------------------------
# 5. Install packages ONCE with install.packages("package_name")
# 6. Load packages EVERY SESSION with library(package_name)
# 7. tidyverse includes many useful packages for data analysis
# 8. dampack provides specialized functions for cost-effectiveness analysis

## 09.03 Data Manipulation -----------------------------------------------------
# 9. Use the pipe operator %>% to chain operations
# 10. dplyr provides intuitive functions: filter(), mutate(), group_by(), summarise()
# 11. Always explore your data before analysis: str(), summary(), View()

## 09.04 Visualization ---------------------------------------------------------
# 12. ggplot2 creates beautiful, publication-ready visualizations
# 13. Basic template: ggplot(data, aes(x, y)) + geom_*() + labs() + theme_*()
# 14. Save plots with ggsave()

## 09.05 Cost-Effectiveness Analysis -------------------------------------------
# 15. Use calculate_icers() from dampack to conduct CEA
# 16. ICERs compare incremental costs to incremental effects
# 17. Dominated strategies are identified automatically (D = strong, ED = weak)
# 18. plot() method creates cost-effectiveness planes
# 19. Always consider all feasible strategies in your analysis

## 09.06 Project Organization --------------------------------------------------
# 20. Follow a consistent folder structure (data/, figures/, analysis/)
# 21. Save processed data with write.csv()
# 22. Document your work with comments


# 10 Next Steps and Resources --------------------------------------------------

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