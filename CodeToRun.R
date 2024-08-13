# ============================================================================ #
#                                 CODE TO RUN                                  #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #
rm(list = ls())
library("dplyr")
library("here")
library("lubridate")
library("pbatR")
library("tableone")
library("ggplot2")
library("stringr")
library("tidyverse")
library("flextable")
library("ftExtra")
library("coloc")
library("readr")
library("tidyverse")

source(here("R_Scripts/Functions.R"))

# Directory where the data is
dir_data    <- 'D:/Projects/VaccineResponse_GWAS/'
dir_ukb     <- 'D:/Projects/VaccineResponse_GWAS/UKBiobank/'
dir_results <- 'D:/Projects/VaccineResponse_GWAS/Results/'

# Load codes
source(here('R_Scripts','1-ImmuneResponse.R'))
source(here('R_Scripts','2-Breakthrough.R'))

# Run the GWAS in the RAP PLATFORM ---------------------------------------------
source(here("R_Scripts","3-ComputePVal.R"))

# Save the results within the "GWAS" folder under the names:
# breakthroughSeverity.txt
# breakthroughSeverity_Validation.txt

# Run FUMA ---------------------------------------------------------------------

# Colocalisation analysis ------------------------------------------------------
source(here("R_Scripts","4-Colocalisation.R"))

# Validation -------------------------------------------------------------------
source(here("R_Scripts","5-Validation.R"))

# SCRIPTS TO MAKE THE TABLES/FIGURES FROM THE PAPER ----------------------------
source(here("R_Scripts","6-CreateTables.R"))

source(here("R_Scripts","7-CreateManhattanPlots.R"))
source(here("R_Scripts","8-CreateValidationPlot.R"))


