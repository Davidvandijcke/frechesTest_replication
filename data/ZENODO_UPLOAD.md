# Instructions for openICPSR Data Archive

This document provides instructions for uploading the replication data to openICPSR.

## Archive Structure

The prepared `data_clean/` folder (48GB) contains:

```
data_clean/
├── README.md                     # Detailed description of contents
├── MANIFEST.txt                  # File inventory and sizes
├── in/
│   ├── sipp_raw/                # 52 SIPP raw zip files
│   ├── sipp_unzipped/           # 10 pre-processed .dta files
│   ├── SIPP_data_dictionaries/  # PDF documentation
│   ├── eora_io_data_2010/       # EORA I-O data by year
│   ├── eora_io_data_2011/
│   ├── ... (through 2021)
│   ├── oecd_io_data/            # OECD I-O .rdata files
│   ├── Eora26Structure.xlsx
│   └── getOECD_IO.R
└── out/
    └── [empty - files created during replication]
```

## Uploading to openICPSR

### 1. Prepare the Archive
```bash
# From the directory containing data_clean/
tar -czf frechet_anova_data.tar.gz data_clean/
# or for better compression:
zip -r frechet_anova_data.zip data_clean/
```

### 2. Create openICPSR Project
1. Log in to https://www.openicpsr.org
2. Click "Share Data" → "Get Started"
3. Select appropriate deposit type
4. Create new project

### 3. Upload Data
1. Upload `frechet_anova_data.tar.gz` (or .zip)
2. Also upload separately:
   - `data_clean/README.md`
   - `data_clean/MANIFEST.txt`

### 4. Add Metadata
- **Title**: "Replication Data for: A Test for Jumps in Metric-Space Conditional Means"
- **Principal Investigator**: Van Dijcke, David (University of Michigan)
- **Description**: 
  ```
  This dataset contains the replication data for "A Test for Jumps in 
  Metric-Space Conditional Means". It includes:
  
  1. SIPP (Survey of Income and Program Participation) data from the U.S. 
     Census Bureau, years 2008-2023, used for occupation/industry transition 
     analysis
  
  2. EORA Global Supply Chain Database input-output tables for years 2010-2021, 
     used for World Bank income classification analysis
  
  3. OECD Input-Output data in pre-processed format
  
  The archive is organized to work directly with the replication code 
  available at: [GitHub URL]
  ```
- **Keywords**: Fréchet test, jump test, regression discontinuity, metric spaces, SIPP, input-output analysis, EORA
- **Geographic Coverage**: Global (for I-O data), United States (for SIPP)
- **Time Period**: 2008-2023 (SIPP), 2010-2021 (EORA)
- **Collection Notes**: Document that this is secondary data compiled from public sources

### 5. Set Access Level
- Choose appropriate access level (likely "Public Use")
- Add any necessary use restrictions

### 6. Link to Code
In the documentation, clearly state:
```
Replication code available at: https://github.com/[your-username]/frechet_anova

Users should:
1. Download and extract this data archive
2. Rename the folder from 'data_clean' to 'data'
3. Place it at the same level as the 'code' directory from GitHub
```

## After Upload

1. **Get DOI**: Once published, openICPSR will assign a DOI
2. **Update Code Repository**:
   - Update `code/README.md` with the openICPSR DOI
   - Update `code/data/README_DATA.md` with the direct link
3. **Test Download**: Verify the archive downloads and extracts correctly

## Citation Format

After receiving DOI:
```
Van Dijcke, David. 2025. "Replication Data for: A Test for Jumps in 
Metric-Space Conditional Means." ICPSR - Interuniversity Consortium 
for Political and Social Research. https://doi.org/10.3886/[DOI]
```