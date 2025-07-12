# Quick Start: Data Setup

## Option 1: Download from openICPSR (Recommended)

1. **Download the data archive** from openICPSR:
   - URL: [To be added after upload]
   - Size: 31GB compressed

2. **Extract and organize**:
   ```bash
   # Extract the archive
   tar -xzf frechet_anova_data.tar.gz
   # or
   unzip frechet_anova_data.zip
   
   # Rename the folder
   mv data_clean data
   ```

3. **Verify directory structure**:
   ```
   your_working_directory/
   ├── code/        # This repository
   ├── data/        # The renamed data folder
   ├── figs/        # Created by code
   └── tabs/        # Created by code
   ```

4. **Run the code**:
   ```r
   setwd("path/to/code")
   source("000_master.R")
   ```

## Option 2: Download from Original Sources

See `data/README_DATA.md` for detailed instructions on downloading from:
- U.S. Census Bureau (SIPP data)
- EORA World MRIO (Input-Output tables)
- OECD.Stat (Alternative I-O data)

Note: This option requires significantly more time and disk space (>200GB).

## Troubleshooting

- **Path errors**: Ensure `data/` is at `../data/` relative to `code/`
- **Memory issues**: The I-O analysis requires ~16GB RAM
- **Missing files**: Check that all subdirectories were extracted properly

For detailed instructions, see `data/README_DATA.md`