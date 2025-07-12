# Data Instructions for Replication

This document provides instructions for obtaining and organizing the data needed for replication.

## Quick Start: Download from openICPSR

**The complete replication data is available at**: https://huggingface.co/datasets/dvdijcke/frechesTest_replication

### Setup Instructions:
1. Download the data archive from huggingface ()
2. Rename it to `data/`
3. Place it at the same level as the `code/` directory (not inside it)

Your directory structure should look like:
```
frechet_anova/
├── code/        # This repository
├── data/        # Data from openICPSR (renamed from data_clean/)
│   ├── in/
│   └── out/
├── figs/        # Will be created by code
└── tabs/        # Will be created by code
```

## Data Contents

The openICPSR archive contains all necessary data files:

### Input Data (`data/in/`)
- **SIPP Data** (2018-2023 only)
  - `sipp_raw/`: 12 raw zip files from Census Bureau (2018-2023)
  - `sipp_unzipped/`: 6 pre-unzipped .dta files (2018-2023)
  - `SIPP_data_dictionaries/`: PDF documentation

- **EORA Input-Output Tables** (2010-2021)
  - `eora_io_data_YYYY/`: Separate folder for each year
  - Each contains trade matrices (T), value added (VA), final demand (FD), etc.

- **OECD Input-Output Data**
  - `oecd_io_data/`: Pre-processed .rdata files

### Output Data (`data/out/`)
- Initially empty (files created during replication)
- `sipp_job_panel.csv` will be created by `processSIPP.R`

## Alternative: Download from Original Sources

If you prefer to download data from original sources:

### 1. SIPP Data
```r
# From the code directory:
source("Applications/10_grabSIPP.R")  # Downloads from Census Bureau
source("Applications/processSIPP.R")   # Processes the data
```

### 2. EORA Input-Output Tables
1. Register at https://worldmrio.com/
2. Download "Eora26" Basic Prices tables for years 2010-2021
3. Extract each year to `data/in/eora_io_data_YYYY/`

### 3. OECD Data
```r
# The script to download OECD data is at:
source("../data/in/getOECD_IO.R")
```

## Important Notes

1. **Path Configuration**: The code expects the data folder to be at `../data/` relative to the code directory. This is set in `000_master.R`.

2. **Memory Requirements**: Processing the full dataset requires ~16GB RAM. For limited memory:
   - Process years sequentially in `WID_PPP.R`
   - Use fewer parallel cores in simulations

3. **Disk Space**: Ensure at least 80GB free space for:
   - 31GB data archive
   - ~20GB for intermediate files during processing
   - Output figures and tables

## Troubleshooting

- **Path errors**: Verify data is in `../data/` relative to code directory
- **Memory errors**: Reduce parallel cores or process data in chunks
- **Missing files**: Check that extraction preserved the directory structure

For questions: dvdijcke@umich.edu