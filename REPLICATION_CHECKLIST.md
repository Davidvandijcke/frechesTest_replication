# Replication Checklist

This checklist helps verify that all analyses in the paper can be successfully replicated.

## Prerequisites
- [ ] R version 4.0 or higher installed
- [ ] All required packages installed (run `source("00_prep.R")`)
- [ ] Data downloaded from openICPSR and organized per `data/README_DATA.md`
  - [ ] Download data archive from: [openICPSR URL - to be added]
  - [ ] Extract and rename `data_clean/` to `data/`
  - [ ] Place `data/` at same level as `code/` directory

## Simulations (Section 3 of paper)

### Table 1: Monte Carlo Simulation Results
- [ ] Run `source("Simulations/frechesTest_simulations.R")`
- [ ] Check output: `tabs/simulation_results_tabular.tex`
- [ ] Verify CSV: `data/out/sims_table_df.csv`

### Figure 1: Power Curves Comparison
- [ ] Run `source("Simulations/frechetANOVA_simulations.R")`
- [ ] Check output: `figs/power_curves_stacked_2x3_N200.png`

## Empirical Applications (Section 4 of paper)

### Application 1: SIPP Occupation Analysis

#### Data Preparation
- [ ] If using openICPSR data: SIPP data is already included in `data/in/`
- [ ] Alternative: Run `source("Applications/10_grabSIPP.R")` to download from Census
- [ ] Run `source("Applications/processSIPP.R")` to process data
- [ ] Verify output: `data/out/sipp_job_panel.csv` exists

#### Analysis
- [ ] Run `source("Applications/occ_anova.R")`

#### Expected Outputs
- [ ] Figure 2: `figs/transition_heatmap.png`
- [ ] Figure 3: `figs/wfh_composition_combined.png`
- [ ] Table 2: `tabs/summary_stats_table.tex`
- [ ] Table 3: `tabs/rdd_balance_tests.tex`

### Application 2: World Bank Income Classification

#### Data Preparation
- [ ] If using openICPSR data: EORA data is already included in `data/in/eora_io_data_YYYY/`
- [ ] Verify EORA data folders exist for years 2010-2021
- [ ] Note: Analysis uses years 2015-2017 by default (configurable in script)

#### Analysis
- [ ] Run `source("Applications/WID_PPP.R")`

#### Expected Outputs
- [ ] Figure 4: Network difference plots
  - [ ] `figs/network_difference_A_matrix.png`
  - [ ] `figs/network_difference_L_matrix.png`
- [ ] Figure 5: RDD plots for scalar outcomes
  - [ ] `figs/rdd_plot_upstream_centrality.png`
  - [ ] `figs/rdd_plot_downstream_centrality.png`
  - [ ] `figs/rdd_plot_manufacturing_intensity.png`
  - [ ] `figs/rdd_plot_services_intensity.png`
- [ ] Table 4: `tabs/scalar_rdd_results_summary.csv`
- [ ] Additional results: `tabs/frechet_test_results_summary.csv`

## Troubleshooting

### Memory Issues
- For World Bank analysis, ensure at least 16GB RAM available
- Process years sequentially if needed by modifying the year range in `WID_PPP.R`

### Missing Packages
- If a package is missing, `00_prep.R` should install it automatically
- For manual installation: `install.packages("package_name")`

### Data Download Issues
- If using openICPSR archive: Ensure complete extraction and correct directory structure
- For manual downloads:
  - SIPP: Check internet connection, Census API may have rate limits
  - EORA: Try downloading during off-peak hours, contact support@worldmrio.com if issues persist

## Final Verification
- [ ] All figures in `figs/` directory match paper
- [ ] All tables in `tabs/` directory match paper
- [ ] No error messages during execution
- [ ] Results are reproducible across multiple runs (set seed is used)

## Notes
- Total runtime: ~4-6 hours depending on system
- Simulations can be run in parallel (default) or serial mode
- Some randomness is expected in bootstrap confidence intervals