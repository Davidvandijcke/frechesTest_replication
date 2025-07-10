# Replication Package: A Test for Jumps in Metric-Space Conditional Means

This repository contains the replication code for "A Test for Jumps in Metric-Space Conditional Means" by David Van Dijcke.

## Overview

The code implements and demonstrates the Fréchet Jump Test, a method for detecting discontinuities in conditional means of random objects valued in metric spaces. The package includes simulations and empirical applications.

## Requirements

### Software
- R version 4.0 or higher
- RStudio (recommended)

### R Packages
The required packages are automatically installed by running `00_prep.R`. Key packages include:
- `frechesTest`: The main package implementing the jump test (included in `frechesTest/` directory)
- `frechet`: For Fréchet regression and ANOVA
- Standard packages: `data.table`, `ggplot2`, `rdrobust`, `rddensity`, `foreach`, `doParallel`

## Repository Structure

```
├── 000_master.R              # Master script - sets paths and parameters
├── 00_prep.R                 # Loads packages and sets themes
├── utils.R                   # Utility functions
│
├── frechesTest/              # R package implementing the Fréchet Jump Test
│   ├── R/                    # Package source code
│   ├── tests/                # Package tests
│   └── man/                  # Package documentation
│
├── Simulations/
│   ├── frechesTest_simulations.R    # Monte Carlo simulations for test performance
│   └── frechetANOVA_simulations.R   # Power curve comparisons with Fréchet ANOVA
│
├── Applications/
│   ├── occ_anova.R           # SIPP occupation/industry transition analysis
│   ├── WID_PPP.R             # World Bank income classification analysis
│   ├── 10_grabSIPP.R         # SIPP data download script
│   └── processSIPP.R         # SIPP data processing
│
├── data/                     # Data directory (see Data section below)
├── figs/                     # Output figures
└── tabs/                     # Output tables
```

## Data Requirements

Due to size constraints, data files are not included in this repository. The code uses publicly available data that can be obtained as follows:

### 1. SIPP Data (for occupation analysis)
- Run `10_grabSIPP.R` to download SIPP data from Census Bureau
- Then run `processSIPP.R` to process the raw data
- Output: `data/out/sipp_job_panel.csv`

### 2. World Bank Data (for income classification analysis)
- Automatically downloaded via WDI package when running `WID_PPP.R`
- EORA Input-Output tables: Download from https://worldmrio.com/
- Place in `data/in/eora_io_data/` with structure:
  ```
  data/in/eora_io_data/
  ├── 2010/
  ├── 2011/
  └── ... (through 2015)
  ```

### 3. Simulation Data
- No external data required - simulations generate synthetic data

## Replication Instructions

### 1. Initial Setup
```r
# Set working directory to the code folder
setwd("path/to/code")

# Run setup
source("000_master.R")
```

### 2. Run Simulations
```r
# Test performance simulations (Table 1 and related figures)
source("frechesTest_simulations.R")

# Power curve comparisons (Figure 1)
source("frechetANOVA_simulations.R")
```

### 3. Run Empirical Applications
```r
# SIPP occupation analysis (Figures 2-3, Tables 2-3)
# First download and process SIPP data:
source("10_grabSIPP.R")
source("processSIPP.R")
# Then run analysis:
source("occ_anova.R")

# World Bank income classification analysis (Figures 4-5, Table 4)
# Ensure EORA data is downloaded first, then:
source("WID_PPP.R")
```

## Output

- **Figures**: All figures are saved to `figs/` directory as PNG/PDF files
- **Tables**: LaTeX tables are saved to `tabs/` directory
- **Results**: Detailed results are saved as RDS files in `figs/raw_results/`

## Notes

1. The `frechesTest` package is included as a subdirectory. To install it:
   ```r
   devtools::install("frechesTest")
   ```

2. Parallel processing is enabled by default in simulations. To run serially (for debugging), set `parallel_sim <- FALSE` in the simulation scripts.

3. Some analyses may take several hours to run, especially the World Bank analysis with multiple years of I-O data.

## Citation

If you use this code, please cite:

```
Van Dijcke, David (2025). "A Test for Jumps in Metric-Space Conditional Means." 
Working Paper, University of Michigan.
```

## Contact

David Van Dijcke  
University of Michigan  
Email: dvdijcke@umich.edu

## License

This code is provided under the MIT License. See LICENSE file for details.