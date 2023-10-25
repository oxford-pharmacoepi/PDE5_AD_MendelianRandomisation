rm(list = ls())
pacman::p_load('readr','devtools','dplyr','TwoSampleMR','LDlinkR','forestplot','tidyr')
devtools::install_github("bar-woolf/TwoStepCisMR")
library(TwoStepCisMR)
pathToData <- 'C:/Users/marta/Desktop/PDE5_AD_MendelianRandomisation/'

# Mendelian randomisation analysis ---------------------------------------------
