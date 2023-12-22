rm(list = ls())
library("magrittr")

pathData <- 'D:/Projects/PDE5_AD_MendelianRandomisation/'
dir.create(paste0(pathData,"Results/"))
pathResults <- paste0(pathData,"Results/")

# Instrument selection
source(here::here("Scripts/InstrumentSelection.R"))

# Mendelian randomisation
dir.create(paste0(pathResults,"MR_Results"))
source(here::here("Scripts/MendelianRandomisation.R"))

# Colocalisation
dir.create(paste0(pathResults,"Colocalisation"))
source(here::here("Scripts/Colocalisation.R"))

# Two-step mendelian randomisation
dir.create(paste0(pathResults,"TwoStepMR"))
source(here::here("Scripts/TwoStepMR.R"))

# Leave-one-out analysis
dir.create(paste0(pathResults,"LeaveOneOutAnalysis"))
source(here::here("Scripts/LeaveOneOutAnalysis.R"))

# Figures and tables from the article
dir.create(paste0(pathResults,"TablesAndFigures"))
source(here::here("Scripts/TablesAndFigures."))