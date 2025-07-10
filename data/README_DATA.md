# Data Instructions for Replication

This document provides detailed instructions for obtaining and organizing the data needed for replication.

## Data Structure

The code expects the following directory structure:

```
data/
├── in/                      # Input data (raw)
│   ├── eora_io_data/        # EORA Input-Output tables
│   │   ├── 2010/
│   │   ├── 2011/
│   │   ├── 2012/
│   │   ├── 2013/
│   │   ├── 2014/
│   │   └── 2015/
│   └── [SIPP files downloaded by 10_grabSIPP.R]
│
├── by/                      # Intermediate data (created by processing scripts)
│
└── out/                     # Processed data
    └── sipp_job_panel.csv   # Created by processSIPP.R
```

## 1. SIPP Data (Survey of Income and Program Participation)

**Source**: U.S. Census Bureau  
**Access**: Public, downloaded via API

### Instructions:
1. Run `10_grabSIPP.R` to automatically download SIPP data files
2. Run `processSIPP.R` to process the raw files
3. This creates `data/out/sipp_job_panel.csv` used in `occ_anova.R`

**Size**: ~500 MB raw, ~100 MB processed

## 2. EORA Input-Output Tables

**Source**: EORA Global Supply Chain Database  
**Website**: https://worldmrio.com/  
**Access**: Free registration required

### Instructions:
1. Register at https://worldmrio.com/
2. Download "Basic prices" tables for years 2010-2015
3. For each year, download the ZIP file containing:
   - The inter-industry transaction matrix (Z matrix)
   - Labels files
4. Extract each year's data into `data/in/eora_io_data/YYYY/`
5. Required files per year:
   - `Eora26_YYYY_bp/Eora26_YYYY_bp_T.txt` (transaction matrix)
   - `Eora26_YYYY_bp/labels_T.txt` (sector/country labels)

**Size**: ~2 GB total (can be processed year-by-year to save space)

### Alternative: Zenodo Data Archive

To facilitate replication without downloading from multiple sources, we provide a Zenodo archive with all necessary data files pre-organized:

**DOI**: [To be added upon publication]  
**Contents**:
- Processed SIPP panel data
- EORA I-O matrices (2010-2015)
- World Bank income thresholds
- All intermediate files

**Instructions for Zenodo data**:
1. Download the archive from Zenodo
2. Extract to the `data/` directory
3. Verify the directory structure matches above

## 3. Other Data Sources

### World Bank Data
- Downloaded automatically via the `WDI` R package
- No manual download needed

### NAICS/SOC Classification Data
- Downloaded automatically from BLS websites in `occ_anova.R`
- No manual download needed

## Data Size Management Strategy

Since GitHub has file size limitations (100 MB per file, 1 GB recommended per repo), we use the following strategy:

1. **Code Repository**: Contains only code and documentation (~10 MB)
2. **Data Archive**: Hosted on Zenodo with DOI for permanent access
3. **Download Scripts**: Provided for users who prefer to download from original sources

## Troubleshooting

### EORA Data Issues
- If the EORA website is slow, try downloading during off-peak hours
- Each year can be processed independently if memory is limited
- Contact EORA support if download links are broken

### SIPP Data Issues  
- The Census API may have rate limits; the script includes appropriate delays
- If downloads fail, check your internet connection and try again
- Some files are large; ensure sufficient disk space

### Memory Issues
- The World Bank analysis with I-O matrices requires ~16GB RAM
- Process years sequentially if memory is limited
- Consider using a high-performance computing cluster for full replication

## Contact for Data Questions

For questions about data access or processing:
- SIPP data: Contact U.S. Census Bureau
- EORA data: Contact support@worldmrio.com
- Replication issues: Contact dvdijcke@umich.edu