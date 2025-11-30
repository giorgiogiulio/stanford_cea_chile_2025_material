


library(lhs) # package for Latin Hypercube Sampling


devtools::install_version("IMIS", version = "0.1", repos = "http://cran.us.r-project.org")
library(IMIS) # package for Incremental Mixture Importance Sampling
library(matrixStats) # package used for summary statistics

# visualization
library(plotrix) # for plotCI function
library(psych) # for pairs.panels function
library(scatterplot3d) # now that we have three inputs to estimate, we'll need higher dimension visualization
library(ggplot2) # general plotting
library(GGally) # ggplot-type correlation plot

# data manipulation
library(dplyr) # for data manipulation

