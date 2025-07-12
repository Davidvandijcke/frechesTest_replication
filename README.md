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
│   ├── frechesTest_simulations.R    # Monte Carlo simulations for test performance -- Table 1 in paper
│   └── frechetANOVA_simulations.R   # Power curve comparisons with Fréchet ANOVA -- Figure 1 in paper
│
├── Applications/
│   ├── occ_anova.R           # SIPP occupation/industry transition analysis -- Figure 2 & 3 in paper
│   ├── WID_PPP.R             # World Bank income classification analysis -- Figure 4 in paper
│   ├── 10_grabSIPP.R         # SIPP data download script
│   └── processSIPP.R         # SIPP data processing
│
├── data/                     # Data directory (see Data section below)
├── figs/                     # Output figures
└── tabs/                     # Output tables
```

## Data Requirements

The complete data required for replication is available from openICPSR at: [URL to be added after upload]

### Data Download and Setup

1. **Download the replication data** from openICPSR
   - The data archive contains a `data_clean/` folder (31GB)
   - Extract this folder and rename it to `data/`
   - Place it in the same directory as this code folder

2. **Verify the data structure**:
   ```
   data/
   ├── in/
   │   ├── sipp_raw/             # Raw SIPP data files (12 zip files, 2018-2023 only)
   │   ├── sipp_unzipped/        # Pre-unzipped SIPP data (6 .dta files, 2018-2023 only)
   │   ├── SIPP_data_dictionaries/
   │   ├── eora_io_data_2010/    # EORA I-O data by year
   │   ├── eora_io_data_2011/
   │   ├── ... (through 2021)
   │   ├── oecd_io_data/         # OECD I-O data files
   │   └── Eora26Structure.xlsx
   └── out/
       └── sipp_job_panel.csv    # Pre-processed SIPP panel (if included)
   ```

### Data Sources (for reference)

The data archive contains processed versions of:

1. **SIPP Data** (Survey of Income and Program Participation)
   - Source: U.S. Census Bureau
   - Years: 2018-2023 (only these years are included and used)
   - Used for occupation/industry transition analysis in Washington state

2. **EORA Input-Output Tables**
   - Source: https://worldmrio.com/
   - Years: 2010-2021
   - Used for World Bank income classification analysis

3. **OECD Input-Output Data**
   - Source: OECD.Stat
   - Pre-processed .rdata files included

### Alternative: Download Raw Data

If you prefer to download the original data sources:

1. **SIPP Data**: Run `Applications/10_grabSIPP.R` (downloads from Census Bureau)
2. **EORA Data**: Register and download from https://worldmrio.com/
3. **OECD Data**: Run the included `data/in/getOECD_IO.R` script

Note: The raw data download will be significantly larger (>200GB) and require additional processing.

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
# First process SIPP data:
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