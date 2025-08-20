#----------------------------------1. Set phenotypes and generate Rdata---------------------------------
rm(list = ls()) 
gc()
####----------------- Set phenotypes ---------###
library(readxl)
code <- read_xlsx("./Database/Disease_Codes_[240112].xlsx", sheet = "EAS")
code <- as.data.frame(code)
# Match the phenotype column in code!!!!
sumstats1_for = c("IBD")
sumstats2_for = c("SLE")
#-------- Set parameters -------#
path_prefix = "./202410/" ### This path must exist
file_name = "IBD_SLE" ### This folder will be created, PLACO results will be placed in the PLACO subfolder
PLACO_path = "./scripts/PLACO/PLACO.RData" ### PLACO package
#------ Run code ------#
source("./scripts/PLACO/1.PLACO_generate_Rdata(no_MHC_removal).R")

###--------------------- PLACO loop --------------------------
ncore = 4 # Number of cores for parallel processing
source("./scripts/PLACO/2.PLACO_calculation_loop.R")