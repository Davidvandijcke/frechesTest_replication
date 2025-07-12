#****************************************************************************************************************************************************

# MASTER SCRIPT: Metric Space Jump Test - TEST VERSION WITH data_clean

# This is a temporary version that uses data_clean instead of data for testing

#****************************************************************************************************************************************************

#### SET OVERALL PARAMETERS ####


#### SET PATHS ####

if (!require("here", character.only=T)) {install.packages("here", dependencies=TRUE)}; require("here")
codeDir <- here::here()
setwd(codeDir) # sets cd to program directory

dir <- sub("/[^/]+$", "", codeDir)# get main directory

# MODIFIED: Use data_clean instead of data
dataIn <- file.path(dir, "data_clean", "in")
dataBy <- file.path(dir, "data_clean", "by")
dataOut <- file.path(dir, "data_clean", "out")
tabs <- file.path(codeDir, "tabs")
figs <- file.path(codeDir, "figs")

# Create output directories if they don't exist
output_dirs <- list(dataBy, dataOut, tabs, figs,
                   file.path(figs, "plots"),
                   file.path(figs, "tables"),
                   file.path(figs, "raw_results"),
                   file.path(figs, "robustness"),
                   file.path(figs, "logs"))

for (d in output_dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

#### LOAD PACKAGES AND SET THEME ####

source(file.path(codeDir, "00_prep.R"))

# Print configuration message
cat("\n")
cat("=====================================\n")
cat("   REPLICATION CODE - TEST MODE\n")
cat("   Using data_clean directory\n")
cat("=====================================\n")
cat("\nPaths configured:\n")
cat("• Code:", codeDir, "\n")
cat("• Data:", dataIn, "\n")
cat("• Output:", dataOut, "\n")
cat("• Figures:", figs, "\n")
cat("• Tables:", tabs, "\n")

# Verify data availability
cat("\nData verification:\n")
sipp_files <- list.files(file.path(dataIn, "sipp_unzipped"), pattern = "\\.dta$")
cat("• SIPP files:", length(sipp_files), "files (2018-2023)\n")

eora_years <- c()
for (y in 2010:2021) {
  if (dir.exists(file.path(dataIn, paste0("eora_io_data_", y)))) {
    eora_years <- c(eora_years, y)
  }
}
cat("• EORA data: Years", paste(range(eora_years), collapse="-"), "\n")

cat("\nReady to run analyses!\n")
cat("=====================================\n\n")


#### RUN ANALYSES ####
source(file.path(codeDir, "utils.R"))


## Process the SIPP data -- takes a moment to finish
source(file.path(codeDir, "Applications", "processSIPP.R"))

# Application 1: WFH
source(file.path(codeDir, "Applications", "occ_anova.R"))

source(file.path(codeDir, "Applications", "WID_PPP.R"))

