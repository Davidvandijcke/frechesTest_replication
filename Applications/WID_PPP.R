# =============================================================================
# World Bank Income Classification RDD Analysis with Fréchet Jump Tests
# Improved Version with Better Structure, Error Handling, and PPP Scaling
# =============================================================================

# --- CONFIGURATION SECTION ---
config <- list(
  # Directories
  data_dir = dataIn, # Placeholder: "path/to/your/input_data_directory"
  output_dir = figs, # Placeholder: "path/to/your/output_figures_tables_directory"
  cache_dir = file.path(figs, "eora_cache"), # Adjusted to use figs as base
  
  # Analysis options
  analysis_mode = "pooled_years", # "single_year" or "pooled_years"
  single_year_choice = 2015, 
  pooled_years = c(2015, 2016, 2017),
  
  max_io_year_available = 2022,
  
  # RDD options
  bandwidth_limit = 20, # In thousands of USD from the cutoff
  min_obs_per_side = 5,
  
  # Fréchet test options
  cv_k_folds = 5,
  cv_n_bw_candidates = 10,
  
  # For run_analysis_with_io_lags
  base_rv_year_for_specific_lags = 2015,
  io_lags_to_test = c(0, 1, 2) # 0 means IO year = RV year
)

# --- SETUP SECTION ---
required_packages <- c(
  "tidyverse", "rdrobust", "rddensity", "WDI", "readxl",
  "countrycode", "fixest", "cli", "rlang", "ggplot2",
  "stringr", "stats", "utils", "grDevices", "scales", "tidyr"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
invisible(lapply(required_packages, install_if_missing))

setup_directories <- function(base_dir, cache_dir_path) {
  subdirs <- c("tables", "plots", "robustness", "raw_results", "data", "logs", "plots/frechet_diff") # Added plots/frechet_diff
  for (subdir in subdirs) {
    dir_path <- file.path(base_dir, subdir)
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  }
  if (!dir.exists(cache_dir_path)) {
    dir.create(cache_dir_path, recursive = TRUE)
    cat("Created cache directory at:", cache_dir_path, "\n")
  }
}

# --- HELPER FUNCTIONS SECTION ---

create_obs_id <- function(iso3c, year_identifier) {
  paste0(iso3c, "_", year_identifier)
}

check_manipulation <- function(data, config_opts) {
  tryCatch(
    {
      if (!"running_var_centered" %in% names(data)) {
        warning("`running_var_centered` not found for manipulation test")
        return(NULL)
      }
      data_to_test <- data[!is.na(data$running_var_centered), ]
      if (nrow(data_to_test) < (config_opts$min_obs_per_side * 4)) {
        warning("Insufficient valid data for manipulation test (rddensity).")
        return(NULL)
      }
      if (sum(data_to_test$running_var_centered < 0, na.rm = TRUE) < 2 ||
          sum(data_to_test$running_var_centered >= 0, na.rm = TRUE) < 2) {
        warning("Not enough obs on both sides of cutoff for rddensity")
        return(NULL)
      }
      rddensity::rddensity(data_to_test$running_var_centered, c = 0)
    },
    error = function(e) {
      warning(paste("Manipulation test (rddensity) failed:", e$message))
      NULL
    }
  )
}

plot_rdd <- function(data, outcome_var, bandwidth = NULL, config_opts) {
  if (!outcome_var %in% names(data)) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = paste("Error: Variable", outcome_var, "not found")) +
        ggplot2::theme_minimal()
    )
  }
  
  plot_df <- data %>%
    dplyr::filter(!is.na(.data[[outcome_var]]), !is.na(running_var_centered))
  
  if (nrow(plot_df) < (config_opts$min_obs_per_side * 2)) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = paste("Insufficient data to plot RDD for", outcome_var)) +
        ggplot2::theme_minimal()
    )
  }
  
  if (is.null(bandwidth)) {
    rd_result_for_plot_bw <- tryCatch(
      rdrobust::rdrobust(plot_df[[outcome_var]], plot_df$running_var_centered, c = 0),
      error = function(e) NULL
    )
    bandwidth <- if (!is.null(rd_result_for_plot_bw) &&
                     !is.na(rd_result_for_plot_bw$bws[1, 1]) &&
                     rd_result_for_plot_bw$bws[1, 1] > 0) {
      rd_result_for_plot_bw$bws[1, 1]
    } else {
      max_abs_rv <- max(abs(plot_df$running_var_centered), na.rm = TRUE)
      default_bw <- if (is.finite(max_abs_rv) && max_abs_rv > 0) max_abs_rv / 3 else 5
      warning(
        paste(
          "Optimal BW for plotting", outcome_var,
          "NA/invalid, using heuristic:", round(default_bw, 2)
        )
      )
      default_bw
    }
  }
  
  plot_window_half_width <- min(
    bandwidth * 1.5,
    max(abs(plot_df$running_var_centered), na.rm = TRUE) * 0.8,
    na.rm = TRUE
  )
  if (!is.finite(plot_window_half_width) || plot_window_half_width <= 0) plot_window_half_width <- 5
  
  plot_data_binned <- plot_df %>%
    dplyr::filter(abs(running_var_centered) <= plot_window_half_width) %>%
    dplyr::mutate(
      bin = cut(running_var_centered,
                breaks = scales::pretty_breaks(n = 20)(c(-plot_window_half_width, plot_window_half_width)),
                include.lowest = TRUE, right = FALSE
      )
    ) %>%
    dplyr::filter(!is.na(bin)) %>%
    dplyr::group_by(bin, treatment) %>%
    dplyr::summarise(
      mean_outcome = mean(!!rlang::sym(outcome_var), na.rm = TRUE),
      se_outcome = stats::sd(!!rlang::sym(outcome_var), na.rm = TRUE) / sqrt(dplyr::n()),
      running_var_mid = mean(running_var_centered, na.rm = TRUE),
      n_obs = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(mean_outcome), n_obs > 0)
  
  if (nrow(plot_data_binned) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = paste("No binned data to plot RDD for", outcome_var)) +
        ggplot2::theme_minimal()
    )
  }
  
  p <- ggplot2::ggplot(plot_data_binned, ggplot2::aes(x = running_var_mid, y = mean_outcome)) +
    ggplot2::geom_point(ggplot2::aes(size = n_obs, color = factor(treatment)), alpha = 0.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_outcome - 1.96 * se_outcome,
        ymax = mean_outcome + 1.96 * se_outcome,
        color = factor(treatment)
      ),
      width = 0.05 * plot_window_half_width, alpha = 0.5
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::geom_smooth(
      data = dplyr::filter(plot_data_binned, treatment == 0 & !is.na(mean_outcome) & n_obs > 1),
      method = "loess", formula = y ~ x, se = TRUE, color = "blue", fill = "blue", alpha = 0.2, na.rm = TRUE
    ) +
    ggplot2::geom_smooth(
      data = dplyr::filter(plot_data_binned, treatment == 1 & !is.na(mean_outcome) & n_obs > 1),
      method = "loess", formula = y ~ x, se = TRUE, color = "darkgreen", fill = "darkgreen", alpha = 0.2, na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = c("0" = "blue", "1" = "darkgreen"), labels = c("Control", "Treated"), name = "Group") +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "N obs in bin") +
    ggplot2::labs(
      title = paste("RDD Plot:", stringr::str_to_title(stringr::str_replace_all(outcome_var, "_", " "))),
      subtitle = paste("Guidance Bandwidth for Bins:", round(bandwidth, 2), "(thousands USD)"),
      x = "GNI pc - Threshold (thousands USD, centered at 0)",
      y = stringr::str_to_title(stringr::str_replace_all(outcome_var, "_", " "))
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"), legend.position = "bottom")
  
  return(p)
}

convert_io_matrices_to_laplacians <- function(matrices_list, method = "weighted", 
                                              laplacian_type = "symmetric") {
  cat("Converting I-O matrices to graph Laplacians (method:", method, 
      ", type:", laplacian_type, ")...\n")
  
  laplacians <- list()
  if (length(matrices_list) == 0) {
    cat("No matrices provided for Laplacian conversion.\n")
    attr(laplacians, "conversion_rate") <- NA
    return(laplacians)
  }
  
  successful_conversions <- 0
  for (obs_id_key in names(matrices_list)) { 
    A <- matrices_list[[obs_id_key]]
    if (is.null(A) || !is.matrix(A) || any(dim(A) == 0) || nrow(A) != ncol(A)) {
      cat("Skipping invalid/non-square matrix for obs_id:", obs_id_key, "\n")
      next
    }
    
    tryCatch({
      W <- A # Using A matrix as the basis for W (weights)
      W[W < 0] <- 0 # Ensure non-negative weights
      if (method == "binary") {
        W <- (W > 1e-6) * 1 # Threshold for binarization
      }
      
      L_val <- NULL
      if (laplacian_type == "out_degree") {
        D_out <- rowSums(W)
        if (method == "normalized" && any(D_out > 0)) { # Random walk: L_rw = I - D_out^(-1) * W
          D_out_inv <- ifelse(D_out > 0, 1/D_out, 0)
          L_val <- diag(nrow(W)) - diag(D_out_inv, nrow=length(D_out_inv)) %*% W
        } else { # Unnormalized: L = D_out - W
          L_val <- diag(D_out, nrow = length(D_out)) - W
        }
      } else if (laplacian_type == "in_degree") {
        D_in <- colSums(W)
        if (method == "normalized" && any(D_in > 0)) { # L_rw_in = I - D_in^(-1) * W^T (less common variant for directed)
          D_in_inv <- ifelse(D_in > 0, 1/D_in, 0)
          L_val <- diag(nrow(W)) - diag(D_in_inv, nrow=length(D_in_inv)) %*% t(W)
        } else { # Unnormalized: L = D_in - W^T
          L_val <- diag(D_in, nrow = length(D_in)) - t(W)
        }
      } else if (laplacian_type == "symmetric") { # Typically L_sym = I - D_total^(-1/2) * (W+W^T)/2 * D_total^(-1/2) or D_total - (W+W^T)/2
        # For A matrix (potentially asymmetric), (W+W^T)/2 makes it symmetric
        W_sym <- (W + t(W)) / 2 
        D_total <- rowSums(W_sym) # For symmetric matrix, rowSums = colSums
        if (method == "normalized" && any(D_total > 0)) { # Symmetric normalized: L_sym = I - D_total^(-1/2) * W_sym * D_total^(-1/2)
          D_sqrt_inv <- ifelse(D_total > 0, 1/sqrt(D_total), 0)
          D_sqrt_inv_mat <- diag(D_sqrt_inv, nrow = length(D_sqrt_inv))
          L_val <- diag(nrow(W_sym)) - D_sqrt_inv_mat %*% W_sym %*% D_sqrt_inv_mat
        } else { # Unnormalized symmetric: L = D_total - W_sym
          L_val <- diag(D_total, nrow = length(D_total)) - W_sym
        }
      } else {
        cat("Unknown laplacian_type:", laplacian_type, "for obs_id:", obs_id_key, ". Skipping.\n")
        next
      }
      
      laplacians[[obs_id_key]] <- L_val
      successful_conversions <- successful_conversions + 1
      
    }, error = function(e) {
      cat("Error converting matrix to Laplacian for obs_id:", obs_id_key, "-", e$message, "\n")
    })
  }
  
  cat("Successfully converted", successful_conversions, "/", length(matrices_list), 
      "matrices to Laplacians\n")
  attr(laplacians, "conversion_rate") <- if (length(matrices_list) > 0) {
    successful_conversions / length(matrices_list)
  } else {
    NA
  }
  
  return(laplacians)
}


# Sector Definitions
EORA_SECTOR_NAMES <- c(
  "Agriculture", "Fishing", "Mining and Quarrying", "Food & Beverages",
  "Textiles and Wearing Apparel", "Wood and Paper",
  "Petroleum, Chemical and Non-Metallic Mineral Products", "Metal Products",
  "Electrical and Machinery", "Transport Equipment", "Other Manufacturing",
  "Recycling", "Electricity, Gas and Water", "Construction",
  "Maintenance and Repair", "Wholesale Trade", "Retail Trade",
  "Hotels and Restraurants", "Transport", "Post and Telecommunications",
  "Financial Intermediation and Business Activities", "Public Administration",
  "Education, Health and Other Services", "Private Households", "Others",
  "Re-export & Re-import"
)
AGG_SECTOR_NAMES <- c("AGR", "MIN", "MAN", "UTL", "CON", "TRD", "FIN_PRO", "PUB_OTH")

map_eora_sector_to_8agg <- function(eora_sector_name) {
  sector_mapping <- list(
    AGR = c("Agriculture", "Fishing"),
    MIN = c("Mining and Quarrying"),
    MAN = c(
      "Food & Beverages", "Textiles and Wearing Apparel", "Wood and Paper",
      "Petroleum, Chemical and Non-Metallic Mineral Products", "Metal Products",
      "Electrical and Machinery", "Transport Equipment", "Other Manufacturing", "Recycling"
    ),
    UTL = c("Electricity, Gas and Water"),
    CON = c("Construction", "Maintenance and Repair"),
    TRD = c(
      "Wholesale Trade", "Retail Trade", "Hotels and Restraurants",
      "Transport", "Post and Telecommunications", "Re-export & Re-import"
    ),
    FIN_PRO = c("Financial Intermediation and Business Activities"),
    PUB_OTH = c(
      "Public Administration", "Education, Health and Other Services",
      "Private Households", "Others"
    )
  )
  for (agg_sector in names(sector_mapping)) {
    if (eora_sector_name %in% sector_mapping[[agg_sector]]) {
      return(agg_sector)
    }
  }
  return(NA_character_)
}
MANUF_AGG_IDX <- which(AGG_SECTOR_NAMES == "MAN")
SERV_AGG_IDX <- which(AGG_SECTOR_NAMES %in% c("TRD", "FIN_PRO", "PUB_OTH"))


# NEW HELPER FUNCTION FOR PLOTTING FRECHET MEAN DIFFERENCES
plot_frechet_mean_difference <- function(left_mean, right_mean, title_suffix = "") {
  if (is.null(left_mean) || is.null(right_mean) ||
      !is.matrix(left_mean) || !is.matrix(right_mean) ||
      any(dim(left_mean) != dim(right_mean)) || all(dim(left_mean) == 0)) {
    cat("Cannot plot Fréchet mean difference: missing, non-matrix, non-conformable, or zero-dimension input for title: ", title_suffix, "\n")
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = paste("Fréchet Mean Difference Plot Error", title_suffix)) +
        ggplot2::theme_minimal()
    )
  }
  
  if (is.null(rownames(left_mean)) || is.null(colnames(left_mean)) ||
      is.null(rownames(right_mean)) || is.null(colnames(right_mean)) ) {
    cat("Error: Matrices for plotting Fréchet mean difference must have dimnames. Title suffix:", title_suffix, "\n")
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = paste("Fréchet Mean Difference Plot Error: Missing Dimnames", title_suffix)) +
        ggplot2::theme_minimal()
    )
  }
  
  diff_matrix <- right_mean - left_mean # Treated - Control
  n_sectors_plot <- nrow(diff_matrix)
  
  plot_df <- tidyr::expand_grid(
    From_Sector_Idx = 1:n_sectors_plot,
    To_Sector_Idx = 1:n_sectors_plot
  ) %>%
    dplyr::mutate(
      Difference = as.vector(diff_matrix), 
      From_Sector = factor(rownames(diff_matrix)[From_Sector_Idx], levels = rownames(diff_matrix)),
      To_Sector   = factor(colnames(diff_matrix)[To_Sector_Idx], levels = colnames(diff_matrix))
    )
  
  max_abs_diff <- max(abs(plot_df$Difference), na.rm = TRUE)
  if (!is.finite(max_abs_diff) || max_abs_diff == 0) max_abs_diff <- 1 
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = From_Sector, y = To_Sector, fill = Difference)) +
    ggplot2::geom_tile(color = "grey50") + 
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      name = "Difference\n(Treated - Control)",
      limits = c(-max_abs_diff, max_abs_diff) 
    ) +
    ggplot2::theme_minimal(base_size = 14) + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 14), 
      axis.text.y = ggplot2::element_text(size = 14),
      # plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 14),
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank()  
    ) +
    ggplot2::coord_fixed() 
  
  return(p)
}


# --- MAIN ANALYSIS FUNCTIONS ---

run_main_analysis <- function(config_opts) {
  cat(
    "=== WORLD BANK INCOME CLASSIFICATION RDD ANALYSIS ===\nAnalysis started:",
    as.character(Sys.time()), "\n"
  )
  setup_directories(config_opts$output_dir, config_opts$cache_dir)
  log_file <- file.path(config_opts$output_dir, "logs",
                        paste0("analysis_log_", format(Sys.Date(), "%Y%m%d"), ".txt")
  )
  sink(log_file, split = TRUE) # Directs console output to log_file and console
  
  tryCatch(
    {
      cat("\n=== STEP 1: Loading World Bank Data ===\n")
      wb_data <- load_world_bank_data() # Now includes PPP indicators
      
      cat("\n=== STEP 2: Preparing Analysis Data (Main Analysis) ===\n")
      # This function will now calculate ppp_scaling_factor
      analysis_data_base <- prepare_analysis_data(wb_data$countries_data, wb_data$thresholds, config_opts)
      if (is.null(analysis_data_base) || nrow(analysis_data_base) == 0) {
        stop("No base analysis data prepared for main analysis. Aborting.")
      }
      
      cat("\n=== STEP 3: Loading EORA I-O Data (Main Analysis) ===\n")
      # process_eora_data will use ppp_scaling_factor from analysis_data_base
      io_data_results <- process_eora_data(analysis_data_base, config_opts)
      if (is.null(io_data_results) || nrow(io_data_results$data) == 0) {
        stop("No I-O data could be processed for main analysis. Aborting.")
      }
      
      cat("\n=== STEP 4: Running Standard RDD (Main Analysis) ===\n")
      rdd_results <- run_scalar_rdd_analysis(io_data_results$data, config_opts)
      
      cat("\n=== STEP 5: Running Fréchet Jump Tests (Main Analysis) ===\n")
      named_running_vars <- stats::setNames(
        io_data_results$data$running_var_centered,
        io_data_results$data$obs_id
      )
      frechet_results <- run_frechet_tests(
        io_data_results$matrices_A, io_data_results$matrices_L,
        named_running_vars, config_opts
      )
      
      cat("\n=== STEP 6: Running Analysis with Fixed RV and I-O Outcome Lags ===\n")
      lagged_io_results <- run_analysis_with_io_lags(
        wb_data, config_opts, # wb_data passed for PPP info
        base_rv_year = config_opts$base_rv_year_for_specific_lags,
        io_outcome_lags = config_opts$io_lags_to_test
      )
      
      cat("\n=== STEP 7: Generating Outputs ===\n")
      generate_all_outputs(
        rdd_results, frechet_results, lagged_io_results,
        io_data_results$data, config_opts
      )
      
      cat("\n=== ANALYSIS COMPLETE ===\nResults saved to:", config_opts$output_dir, "\n")
    },
    error = function(e) {
      cat("\n!!! ERROR IN MAIN ANALYSIS !!!\nError message:", e$message, "\nTraceback:\n")
      print(sys.calls())
    },
    finally = {
      sink() # Reset output to console only
    }
  )
}
load_world_bank_data <- function() {
  wb_thresholds <- data.frame(
    year = 2010:2024,
    low_income_threshold = c(1005, 1025, 1035, 1045, 1045, 1025, 1005, 995, 995, 1025, 1035, 1045, 1085, 1135, 1145),
    high_income_threshold = c(12275, 12475, 12615, 12745, 12735, 12475, 12236, 12056, 12376, 12536, 12696, 13205, 13845, 14005, 14005)
  )
  
  # Define all WDI indicators that are intended to be fetched
  wdi_indicators_to_fetch <- c(
    "NY.GNP.PCAP.CD", "NY.GDP.PCAP.CD", "SP.POP.TOTL", "NE.TRD.GNFS.ZS", # Original
    "PA.NUS.PPP",        # PPP conversion factor, GDP (LCU per international $)
    "PA.NUS.FCRF",       # Official exchange rate (LCU per US$, period average)
    "NY.GDP.MKTP.PP.CD", # GDP, PPP (current international $)
    "NY.GDP.MKTP.PP.KD"  # GDP, PPP (constant 2017 international $)
  )
  
  countries_data_raw <- WDI::WDI(
    indicator = wdi_indicators_to_fetch,
    start = min(wb_thresholds$year) - 5,
    end = max(wb_thresholds$year),
    extra = TRUE
  )
  
  # Check if WDI call returned any data at all
  if (is.null(countries_data_raw) || nrow(countries_data_raw) == 0) {
    stop("WDI call returned no data. Check network connection, WDI query parameters, or WDI service status.")
  }
  
  # Ensure all requested WDI indicators are present as columns in the raw fetched data.
  # If a WDI indicator was not returned (e.g., not available), add it as a column of NAs.
  # This prevents dplyr::rename from failing if an "old" name is missing.
  for (indicator_code in wdi_indicators_to_fetch) {
    if (!indicator_code %in% names(countries_data_raw)) {
      warning(paste("WDI indicator '", indicator_code, "' was not found in the data fetched from World Bank.",
                    " It will be added as a column of NAs. Subsequent calculations involving this indicator will likely result in NAs or default values.", sep=""))
      countries_data_raw[[indicator_code]] <- NA_real_
    }
  }
  
  # Perform renames and initial filtering
  countries_data_processed <- countries_data_raw %>%
    # It's crucial that NY.GNP.PCAP.CD exists for this filter.
    # The loop above ensures it does, even if as NAs.
    dplyr::filter(!is.na(NY.GNP.PCAP.CD)) %>% 
    dplyr::rename(
      gni_per_capita = NY.GNP.PCAP.CD,
      gdp_per_capita = NY.GDP.PCAP.CD,
      population = SP.POP.TOTL,
      trade_gdp = NE.TRD.GNFS.ZS,
      ppp_lcu_per_intl_dollar = PA.NUS.PPP,
      exch_lcu_per_usd = PA.NUS.FCRF,
      gdp_ppp_curr_intl = NY.GDP.MKTP.PP.CD,
      gdp_ppp_const_2017_intl = NY.GDP.MKTP.PP.KD
    ) %>%
    dplyr::mutate(iso3c = countrycode::countrycode(country, "country.name", "iso3c")) %>%
    dplyr::filter(!is.na(iso3c), nchar(iso3c) == 3) # Filter out non-countries/unmatched
  
  # After processing, one final check for the specifically problematic columns
  # to ensure they are indeed present before returning.
  if (!"ppp_lcu_per_intl_dollar" %in% names(countries_data_processed)) {
    stop("FATAL ERROR: Column 'ppp_lcu_per_intl_dollar' is unexpectedly missing after WDI data processing. Original WDI code was 'PA.NUS.PPP'.")
  }
  if (!"exch_lcu_per_usd" %in% names(countries_data_processed)) {
    stop("FATAL ERROR: Column 'exch_lcu_per_usd' is unexpectedly missing after WDI data processing. Original WDI code was 'PA.NUS.FCRF'.")
  }
  
  return(list(countries_data = countries_data_processed, thresholds = wb_thresholds))
}
prepare_analysis_data <- function(countries_data, wb_thresholds, config_opts) {
  cat("Mode:", config_opts$analysis_mode, "\n")
  processed_data_list <- list()
  
  target_rv_years <- if (config_opts$analysis_mode == "single_year") {
    config_opts$single_year_choice
  } else {
    config_opts$pooled_years
  }
  
  for (rv_year_val in target_rv_years) {
    io_year_val <- if (config_opts$analysis_mode == "pooled_years") rv_year_val + 2 else rv_year_val
    cat("  Preparing for RV year:", rv_year_val, "-> Target I-O year for this instance:", io_year_val, "\n")
    
    if (io_year_val > config_opts$max_io_year_available && config_opts$analysis_mode == "pooled_years") {
      cat("    Skipping RV year", rv_year_val, "for pooled: Resulting I-O year", io_year_val,
          "> max available", config_opts$max_io_year_available, "\n")
      next
    }
    
    # df_loop gets GNI etc. for rv_year_val. It also gets PPP data for rv_year_val from countries_data.
    df_loop <- countries_data %>%
      dplyr::filter(year == rv_year_val) %>% 
      dplyr::left_join(wb_thresholds %>% dplyr::filter(year == rv_year_val), by = "year") %>%
      dplyr::mutate(
        io_target_year = io_year_val, 
        running_var_source_year = rv_year_val
      )
    
    if (nrow(df_loop) > 0) {
      processed_data_list[[as.character(rv_year_val)]] <- df_loop
    } else {
      cat("    No country data for RV year:", rv_year_val, "\n")
    }
  }
  
  if (length(processed_data_list) == 0) {
    warning("No data prepared in prepare_analysis_data initial loop.")
    return(dplyr::tibble())
  }
  
  analysis_data_base <- dplyr::bind_rows(processed_data_list)
  
  # ppp_subset_for_join contains PPP data for ALL years from the original countries_data
  ppp_subset_for_join <- countries_data %>%
    dplyr::select(
      iso3c, year, 
      ppp_lcu_per_intl_dollar, exch_lcu_per_usd,
      gdp_ppp_curr_intl, gdp_ppp_const_2017_intl
    ) %>%
    dplyr::rename(data_year_for_ppp = year) # This 'year' is the year of the PPP data
  
  # Join analysis_data_base (GNI for rv_year) with ppp_subset_for_join (PPP for io_target_year)
  # This will create .x and .y columns for the PPP variables because they exist in both
  analysis_data_with_ppp_raw_suffixed <- analysis_data_base %>%
    dplyr::left_join(
      ppp_subset_for_join,
      by = c("iso3c", "io_target_year" = "data_year_for_ppp")
    )
  
  # --- FIX: Rename .y columns (from ppp_subset_for_join, for io_target_year) to desired names ---
  # --- and remove the .x columns (from analysis_data_base, for rv_year_val) ---
  analysis_data_with_ppp_renamed <- analysis_data_with_ppp_raw_suffixed %>%
    dplyr::select(
      -dplyr::any_of(c("ppp_lcu_per_intl_dollar.x", "exch_lcu_per_usd.x", 
                       "gdp_ppp_curr_intl.x", "gdp_ppp_const_2017_intl.x")) # Remove .x cols
    ) %>%
    dplyr::rename(
      ppp_lcu_per_intl_dollar = ppp_lcu_per_intl_dollar.y, # Keep .y and rename
      exch_lcu_per_usd = exch_lcu_per_usd.y,
      gdp_ppp_curr_intl = gdp_ppp_curr_intl.y,
      gdp_ppp_const_2017_intl = gdp_ppp_const_2017_intl.y
    )
  
  # --- DEBUG 3 (REVISED CHECK) ---
  cat("\n--- DEBUG POINT 3 REVISED (prepare_analysis_data): after renaming .y columns ---\n")
  cat("Names in analysis_data_with_ppp_renamed:\n")
  print(names(analysis_data_with_ppp_renamed))
  cat("Structure of analysis_data_with_ppp_renamed (relevant columns):\n")
  if ("ppp_lcu_per_intl_dollar" %in% names(analysis_data_with_ppp_renamed) && "exch_lcu_per_usd" %in% names(analysis_data_with_ppp_renamed)) {
    print(str(analysis_data_with_ppp_renamed %>% 
                dplyr::select(iso3c, year, io_target_year, 
                              ppp_lcu_per_intl_dollar, exch_lcu_per_usd), list.len=5))
    cat("Type of ppp_lcu_per_intl_dollar:", class(analysis_data_with_ppp_renamed$ppp_lcu_per_intl_dollar), "\n")
    cat("Type of exch_lcu_per_usd:", class(analysis_data_with_ppp_renamed$exch_lcu_per_usd), "\n")
  } else {
    cat("CRITICAL AFTER RENAME: 'ppp_lcu_per_intl_dollar' or 'exch_lcu_per_usd' STILL NOT FOUND correctly.\n")
  }
  cat("--- END DEBUG POINT 3 REVISED ---\n\n")
  
  
  # Now, the mutate should work with the correctly named columns for io_target_year
  analysis_data_interim <- analysis_data_with_ppp_renamed %>%
    dplyr::mutate(
      temp_pli_usd_per_intl_dollar = as.numeric(ppp_lcu_per_intl_dollar) / as.numeric(exch_lcu_per_usd),
      temp_deflator_intl_dollar_base2017 = as.numeric(gdp_ppp_curr_intl) / as.numeric(gdp_ppp_const_2017_intl),
      initial_ppp_scaling_factor = (1 / temp_pli_usd_per_intl_dollar) / temp_deflator_intl_dollar_base2017
    )
  
  num_defaulted_scaling_factor <- sum(is.na(analysis_data_interim$initial_ppp_scaling_factor) |
                                        !is.finite(analysis_data_interim$initial_ppp_scaling_factor) |
                                        analysis_data_interim$initial_ppp_scaling_factor <= 0, na.rm = TRUE)
  
  analysis_data_final <- analysis_data_interim %>%
    dplyr::mutate(
      ppp_scaling_factor = ifelse(is.na(initial_ppp_scaling_factor) | !is.finite(initial_ppp_scaling_factor) | initial_ppp_scaling_factor <= 0,
                                  1.0, initial_ppp_scaling_factor),
      running_var = gni_per_capita - high_income_threshold, 
      treatment = as.numeric(running_var >= 0),
      running_var_centered = running_var / 1000,
      analysis_year = running_var_source_year 
    ) %>%
    dplyr::filter(
      !is.na(running_var),
      !is.na(low_income_threshold), 
      abs(running_var_centered) <= config_opts$bandwidth_limit,
      !stringr::str_detect(country, stringr::regex("Euro area|European Union|World|income|OECD|Africa|Asia|Pacific|Latin America|Caribbean|Europe|IDA|IBRD|HIPC|dividend", ignore_case = TRUE))
    ) %>%
    dplyr::mutate(obs_id = create_obs_id(iso3c, io_target_year)) %>%
    dplyr::distinct(obs_id, .keep_all = TRUE) %>%
    dplyr::select(-dplyr::starts_with("temp_"), -initial_ppp_scaling_factor) 
  
  cat("  PPP scaling: ", num_defaulted_scaling_factor, " observations (out of ", nrow(analysis_data_interim), 
      ") will use a default scaling factor of 1.0 due to missing/invalid PPP data or calculation issues for their respective io_target_year.\n")
  cat("  NOTE: Applying PPP scaling to underlying monetary flows (Z and X matrices from EORA). The technical coefficient matrix 'A' (and derived measures like multipliers or Fréchet tests on A) will likely be invariant to this scaling, as A_ij = (factor*Z_ij) / (factor*X_j) = Z_ij / X_j. This scaling aims to make the magnitudes of monetary flows comparable in constant 2017 international dollars.\n")
  
  cat("Prepared base data in prepare_analysis_data:", nrow(analysis_data_final), "obs.\n")
  if (nrow(analysis_data_final) > 0) {
    cat("  RV yrs (running_var_source_year):", paste(sort(unique(analysis_data_final$analysis_year)), collapse = ", "), "\n")
    cat("  IO yrs (io_target_year, for I-O data & PPP):", paste(sort(unique(analysis_data_final$io_target_year)), collapse = ", "), "\n")
    cat("  Sample of PPP scaling factors (first 5 non-default where factor != 1.0):", 
        paste(round(head(analysis_data_final$ppp_scaling_factor[analysis_data_final$ppp_scaling_factor != 1.0 & !is.na(analysis_data_final$ppp_scaling_factor)], 5), 4), collapse=", "), "\n")
  }
  return(analysis_data_final)
}


process_eora_data <- function(analysis_data_with_ids_and_target_io, config_opts) {
  if (nrow(analysis_data_with_ids_and_target_io %||% data.frame()) == 0) {
    cat("No analysis data provided to process_eora_data.\n")
    return(NULL)
  }
  # analysis_data_with_ids_and_target_io now includes ppp_scaling_factor
  
  unique_io_outcome_years <- sort(unique(analysis_data_with_ids_and_target_io$io_target_year))
  cat(
    "Unique I-O outcome years to process based on input data:",
    paste(unique_io_outcome_years, collapse = ", "), "\n"
  )
  
  all_measures_dfs <- list()
  all_matrices_A <- list()
  all_matrices_L <- list()
  
  for (current_io_year_val in unique_io_outcome_years) {
    cat("\nProcessing EORA data for I-O outcome year:", current_io_year_val, "\n")
    current_eora_data_for_year <- tryCatch(
      load_eora_data_for_year(config_opts$data_dir, current_io_year_val),
      error = function(e) {
        cat("FATAL: Failed to load EORA data for year", current_io_year_val, ":", e$message, "\n")
        NULL
      }
    )
    if (is.null(current_eora_data_for_year)) {
      cat("Skipping all processing for I-O year", current_io_year_val, "due to EORA loading error.\n")
      next
    }
    
    analysis_subset_for_this_io_year <- analysis_data_with_ids_and_target_io %>%
      dplyr::filter(io_target_year == current_io_year_val)
    
    if (nrow(analysis_subset_for_this_io_year) == 0) {
      cat("No observations in the analysis plan target I-O year", current_io_year_val, "\n")
      next
    }
    
    # analysis_subset_for_this_io_year contains ppp_scaling_factor
    processed_country_data_for_year <- process_io_for_one_outcome_year(
      analysis_subset_for_this_io_year, # This now includes ppp_scaling_factor
      current_eora_data_for_year,
      config_opts
    )
    
    if (!is.null(processed_country_data_for_year$measures_df) &&
        nrow(processed_country_data_for_year$measures_df) > 0) {
      all_measures_dfs[[as.character(current_io_year_val)]] <- processed_country_data_for_year$measures_df
    }
    all_matrices_A <- c(all_matrices_A, processed_country_data_for_year$matrices_A)
    all_matrices_L <- c(all_matrices_L, processed_country_data_for_year$matrices_L)
  }
  
  if (length(all_measures_dfs) == 0) {
    cat("No I-O measures could be computed for any year after EORA processing.\n")
    return(NULL)
  }
  final_io_measures_df_all_years <- dplyr::bind_rows(all_measures_dfs)
  
  # Need to ensure ppp_scaling_factor is also in final_combined_data_with_io if used later
  # It's part of analysis_data_with_ids_and_target_io, so it should be carried by the join
  # if obs_id is the only key.
  final_combined_data_with_io <- dplyr::inner_join(
    analysis_data_with_ids_and_target_io, # This has ppp_scaling_factor
    final_io_measures_df_all_years, # This has IO measures
    by = "obs_id" # Ensure this is the only join key if there are duplicated columns
    # or select columns from final_io_measures_df_all_years carefully
  )
  
  cat(
    "Total observations with successfully processed and merged I-O data:",
    nrow(final_combined_data_with_io), "\n"
  )
  if (nrow(final_combined_data_with_io) == 0) {
    cat("No observations remain after merging I-O measures back to analysis plan.\n")
    return(NULL)
  }
  
  valid_obs_ids_in_final_df <- final_combined_data_with_io$obs_id
  all_matrices_A_aligned <- all_matrices_A[intersect(names(all_matrices_A), valid_obs_ids_in_final_df)]
  all_matrices_L_aligned <- all_matrices_L[intersect(names(all_matrices_L), valid_obs_ids_in_final_df)]
  
  ids_in_A_aligned <- names(all_matrices_A_aligned)
  ids_in_L_aligned <- names(all_matrices_L_aligned)
  common_matrix_ids_after_align <- intersect(ids_in_A_aligned, ids_in_L_aligned)
  
  final_combined_data_really_final <- final_combined_data_with_io %>%
    dplyr::filter(obs_id %in% common_matrix_ids_after_align)
  
  all_matrices_A_final <- all_matrices_A_aligned[common_matrix_ids_after_align]
  all_matrices_L_final <- all_matrices_L_aligned[common_matrix_ids_after_align]
  
  cat(
    "Final, strictly aligned data: ", nrow(final_combined_data_really_final), " obs with ",
    length(all_matrices_A_final), " A-matrices and ", length(all_matrices_L_final), " L-matrices.\n"
  )
  
  return(list(
    data = final_combined_data_really_final, # This data frame contains ppp_scaling_factor
    matrices_A = all_matrices_A_final,
    matrices_L = all_matrices_L_final
  ))
}

load_eora_data_for_year <- function(base_data_input_dir, year) {
  data_dir_for_specific_year <- file.path(base_data_input_dir, paste0("eora_io_data_", year))
  cat(
    "Attempting to load EORA data for year", year, "from directory:",
    data_dir_for_specific_year, "\n"
  )
  if (!dir.exists(data_dir_for_specific_year)) {
    stop(paste("EORA data directory not found for year", year, ":", data_dir_for_specific_year))
  }
  all_files_in_year_dir <- list.files(
    data_dir_for_specific_year,
    pattern = "\\.txt$",
    full.names = TRUE, ignore.case = TRUE
  )
  if (length(all_files_in_year_dir) == 0) {
    stop(paste("No .txt files in EORA dir for year", year, ":", data_dir_for_specific_year))
  }
  
  eora_data <- list(year = year)
  prefix_pattern <- paste0("Eora26_", year, "_bp_") # Assuming Eora26 naming convention
  matrix_suffices <- list(T_matrix = "T\\.txt$", Q_matrix = "Q\\.txt$")
  
  for (mat_name in names(matrix_suffices)) {
    # Try exact prefix first, then allow for case variations or slight deviations if needed
    full_pattern_regex <- paste0(tolower(prefix_pattern), tolower(matrix_suffices[[mat_name]]))
    
    # Search for files matching the pattern, ignoring case for basename
    candidate_files <- all_files_in_year_dir[
      grepl(full_pattern_regex, tolower(basename(all_files_in_year_dir)))
    ]
    
    # If not found, try with original case (as some systems are case-sensitive)
    if(length(candidate_files) == 0) {
      full_pattern_regex_orig_case <- paste0(prefix_pattern, matrix_suffices[[mat_name]])
      candidate_files <- all_files_in_year_dir[
        grepl(full_pattern_regex_orig_case, basename(all_files_in_year_dir))
      ]
    }
    
    
    file_path_to_load <- NA_character_
    if (length(candidate_files) == 1) {
      file_path_to_load <- candidate_files[1]
    } else if (length(candidate_files) > 1) {
      warning(
        paste(
          "Multiple files for", mat_name, "yr", year, "pattern", full_pattern_regex, ":",
          paste(basename(candidate_files), collapse = ", "), ". Using 1st"
        )
      )
      file_path_to_load <- candidate_files[1]
    }
    
    if (!is.na(file_path_to_load) && file.exists(file_path_to_load)) {
      eora_data[[mat_name]] <- utils::read.table(
        file_path_to_load,
        sep = "\t", header = FALSE,
        comment.char = "", quote = "" # Important for EORA files
      )
      cat("Loaded", mat_name, "for", year, "from", basename(file_path_to_load), "\n")
    } else {
      if (mat_name == "T_matrix") { # T_matrix is essential
        stop(
          paste(
            mat_name, "not found for yr", year, ".Pattern '", prefix_pattern,
            matrix_suffices[[mat_name]], "' (or variants) in:", data_dir_for_specific_year
          )
        )
      }
      # Q_matrix might be estimated later if not found
      warning(
        paste(
          mat_name, "not found for yr", year, ".Pattern '", prefix_pattern,
          matrix_suffices[[mat_name]], "' (or variants). May estimate if Q_matrix."
        )
      )
      eora_data[[mat_name]] <- NULL
    }
  }
  return(create_eora_labels(eora_data, data_dir_for_specific_year, year))
}

create_eora_labels <- function(eora_data, data_dir_for_specific_year, year_val) {
  if (is.null(eora_data$T_matrix)) stop("T_matrix is missing for labels, yr ", year_val)
  actual_rows_T_matrix <- nrow(eora_data$T_matrix)
  cat("  Label Creation Yr:", year_val, "\n    T_matrix rows:", actual_rows_T_matrix, "\n")
  
  # Try year-specific label file first, then generic
  label_file_pattern_year_specific <- paste0("labels_T_", year_val, "\\.txt$")
  label_file_pattern_generic <- "labels_T\\.txt$" # Generic name used by some EORA versions
  
  label_files <- list.files(
    data_dir_for_specific_year,
    pattern = label_file_pattern_year_specific,
    full.names = TRUE, ignore.case = TRUE
  )
  if (length(label_files) == 0) { # If year-specific not found, try generic
    label_files <- list.files(
      data_dir_for_specific_year,
      pattern = label_file_pattern_generic,
      full.names = TRUE, ignore.case = TRUE
    )
  }
  
  if (length(label_files) > 0) {
    label_file_to_use <- label_files[1] # Use the first match
    if (length(label_files) > 1) {
      cat("Multiple label files found (e.g., year-specific and generic), using:", basename(label_file_to_use), "\n")
    }
    tryCatch(
      {
        labels_df <- utils::read.table(
          label_file_to_use,
          sep = "\t", header = FALSE,
          stringsAsFactors = FALSE, fill = TRUE, quote = "", comment.char = ""
        )
        cat(
          "    Label file", basename(label_file_to_use), "loaded:",
          nrow(labels_df), "rows &", ncol(labels_df), "cols.\n"
        )
        # Basic validation of label file structure
        if (nrow(labels_df) != actual_rows_T_matrix) {
          stop(
            paste(
              "MISMATCH: T_matrix rows (", actual_rows_T_matrix,
              ") != label file rows (", nrow(labels_df), ") for year", year_val,
              "using file", basename(label_file_to_use)
            )
          )
        }
        if (ncol(labels_df) < 4) { # EORA labels usually have at least country, sector name, code, description
          stop(paste("Label file", basename(label_file_to_use), "for year", year_val, "has fewer than 4 columns."))
        }
        # Assuming EORA standard: Col 1 is country ISO3C, Col 4 is sector name/description
        eora_data$region_labels <- data.frame(region = as.character(labels_df[, 1]), stringsAsFactors = FALSE)
        eora_data$sector_labels <- data.frame(sector = as.character(labels_df[, 4]), stringsAsFactors = FALSE) # Or V3 if that's sector name
        cat("    Found", length(unique(eora_data$region_labels$region)), "unique regions in labels.\n")
        cat("Successfully created labels for year", year_val, "from", basename(label_file_to_use), "\n")
        return(eora_data)
      },
      error = function(e) {
        stop(paste(
          "Error processing label file", basename(label_file_to_use),
          "for year", year_val, ":", e$message
        ))
      }
    )
  } else {
    stop(
      paste(
        "CRITICAL: Label file (labels_T_YYYY.txt or labels_T.txt) not found in",
        data_dir_for_specific_year, "for year", year_val
      )
    )
  }
}

process_io_for_one_outcome_year <- function(analysis_subset, eora_data_for_this_year, config_opts_for_cache) {
  # analysis_subset now contains ppp_scaling_factor
  n_obs_for_year <- nrow(analysis_subset)
  io_matrices_A_year <- list()
  io_matrices_L_year <- list()
  measures_list_year <- list()
  successful_count_year <- 0
  
  cache_path_year <- file.path(config_opts_for_cache$cache_dir,
                               paste0("processed_io_", eora_data_for_this_year$year)
  )
  if (!dir.exists(cache_path_year)) dir.create(cache_path_year, recursive = TRUE)
  
  cat("Processing", n_obs_for_year, "observations for I-O year", eora_data_for_this_year$year, "\n")
  
  for (i in 1:n_obs_for_year) {
    current_obs_id <- analysis_subset$obs_id[i]
    country_iso_code <- analysis_subset$iso3c[i]
    current_ppp_scaling_factor <- analysis_subset$ppp_scaling_factor[i] # Get the scaling factor
    
    cache_file_name <- paste0(country_iso_code, "_io_", eora_data_for_this_year$year, ".Rds")
    cache_file_full_path <- file.path(cache_path_year, cache_file_name)
    
    if (i %% max(1, round(n_obs_for_year / 10, 0)) == 0 || i == n_obs_for_year || i == 1) {
      cat(
        "  Processing obs", i, "/", n_obs_for_year,
        "(id:", current_obs_id, ", country:", country_iso_code, 
        ", PPP factor:", round(current_ppp_scaling_factor,4) ,")\n"
      )
    }
    
    if (file.exists(cache_file_full_path)) {
      # Potentially add a check here: if ppp_scaling_factor logic changed, invalidate cache or re-tag.
      # For now, assume cache is valid if it exists.
      cat("    Loading cached I-O for", country_iso_code, "year", eora_data_for_this_year$year, "\n")
      cached_result <- readRDS(cache_file_full_path)
      # Simple check, could be enhanced with a version/parameter hash
      if (!is.null(cached_result$A_matrix) && !is.null(cached_result$measures_row) && !is.null(cached_result$L_matrix)) {
        io_matrices_A_year[[current_obs_id]] <- cached_result$A_matrix
        io_matrices_L_year[[current_obs_id]] <- cached_result$L_matrix # Retrieve L_matrix
        measures_for_this_obs <- cached_result$measures_row
        measures_for_this_obs$obs_id <- current_obs_id # Ensure obs_id is set correctly
        measures_list_year[[current_obs_id]] <- measures_for_this_obs
        successful_count_year <- successful_count_year + 1
        next
      } else {
        warning("Cached file for ", country_iso_code, " year ",
                eora_data_for_this_year$year, " invalid or incomplete. Reprocessing."
        )
      }
    }
    
    # Pass ppp_scaling_factor to extract_country_io_data
    country_io_extraction_result <- extract_country_io_data(
      eora_data_for_this_year, 
      country_iso_code,
      current_ppp_scaling_factor # Pass the factor
    )
    
    if (!is.null(country_io_extraction_result$A_matrix)) {
      io_matrices_A_year[[current_obs_id]] <- country_io_extraction_result$A_matrix
      
      # Compute measures using the (potentially PPP-scaled) A_matrix
      # Note: A_matrix itself is likely invariant to PPP scaling of Z and X.
      io_measures_computed <- compute_io_measures_revised(
        country_io_extraction_result$A_matrix, MANUF_AGG_IDX, SERV_AGG_IDX
      )
      io_matrices_L_year[[current_obs_id]] <- io_measures_computed$L_matrix # This is L_inv
      
      measures_to_store_for_df <- io_measures_computed[
        setdiff(names(io_measures_computed), c("A_matrix", "L_matrix")) # Exclude matrices from df
      ]
      
      current_measures_row_for_df <- data.frame(
        obs_id = current_obs_id, as.list(measures_to_store_for_df), stringsAsFactors = FALSE
      )
      measures_list_year[[current_obs_id]] <- current_measures_row_for_df
      
      # Cache the results
      measures_for_cache <- data.frame(as.list(measures_to_store_for_df), stringsAsFactors = FALSE)
      cache_data_to_save <- list(
        A_matrix = country_io_extraction_result$A_matrix,
        L_matrix = io_measures_computed$L_matrix, # Cache L_matrix (L_inv)
        measures_row = measures_for_cache
        # Optionally, store ppp_scaling_factor used, for cache validation later
        # ppp_factor_used = current_ppp_scaling_factor 
      )
      saveRDS(cache_data_to_save, file = cache_file_full_path)
      cat(
        "    Saved processed I-O for", country_iso_code, "year",
        eora_data_for_this_year$year, "to cache.\n"
      )
      successful_count_year <- successful_count_year + 1
    } else {
      cat(
        "    Failed I-O extraction for obs_id:", current_obs_id,
        "(Country:", country_iso_code, ") Status:",
        country_io_extraction_result$status, "\n"
      )
    }
  }
  cat(
    "Successfully processed I-O for", successful_count_year, "/", n_obs_for_year,
    "obs for I-O year", eora_data_for_this_year$year, "\n"
  )
  
  final_measures_df_year <- if (length(measures_list_year) > 0) {
    dplyr::bind_rows(measures_list_year)
  } else {
    # Create an empty tibble with expected column names if no successful processing
    # This ensures bind_rows in the calling function doesn't fail with empty lists
    # Get names from a successful computation or define them explicitly
    # For simplicity, returning an empty tibble if measures_list_year is empty.
    # The calling function handles empty all_measures_dfs.
    dplyr::tibble() 
  }
  
  return(list(
    measures_df = final_measures_df_year,
    matrices_A = io_matrices_A_year,
    matrices_L = io_matrices_L_year # Return L_matrices (L_inv)
  ))
}


extract_country_io_data <- function(eora_data, country_iso3c, ppp_scaling_factor) {
  if (is.null(eora_data$region_labels) || is.null(eora_data$sector_labels)) {
    return(list(A_matrix = NULL, status = paste("Region/sector labels missing for EORA yr", eora_data$year)))
  }
  
  country_indices_in_labels <- which(eora_data$region_labels$region == country_iso3c)
  if (length(country_indices_in_labels) == 0) {
    return(list(A_matrix = NULL, status = paste("Country", country_iso3c, "not found in EORA labels for yr", eora_data$year)))
  }
  
  country_start_idx <- min(country_indices_in_labels)
  n_sectors_for_this_country <- length(country_indices_in_labels)
  country_end_idx <- country_start_idx + n_sectors_for_this_country - 1
  cat("    Extracting for", country_iso3c, ": Label idx", country_start_idx, "-", country_end_idx, 
      "(", n_sectors_for_this_country, "sectors labels), PPP factor:", round(ppp_scaling_factor,4), "\n")
  
  if (country_end_idx > nrow(eora_data$T_matrix) || country_end_idx > ncol(eora_data$T_matrix)) {
    return(list(A_matrix = NULL, status = paste("Indices for", country_iso3c, "out of bounds for T_matrix.")))
  }
  
  # Z_country_native is the block of T_matrix for domestic transactions
  Z_country_native <- as.matrix(eora_data$T_matrix)[country_start_idx:country_end_idx, country_start_idx:country_end_idx]
  X_country_native <- NULL # Total output vector
  
  # Attempt to get X_country_native from Q_matrix (total sector output)
  if (!is.null(eora_data$Q_matrix) && ncol(eora_data$Q_matrix) > 0 && nrow(eora_data$Q_matrix) > 0) {
    Q_matrix_full <- as.matrix(eora_data$Q_matrix)
    # Q_matrix can be a row vector or column vector. Check dimensions.
    if (nrow(Q_matrix_full) == 1 && ncol(Q_matrix_full) >= country_end_idx) { # Row vector
      X_country_native <- as.numeric(Q_matrix_full[1, country_start_idx:country_end_idx])
    } else if (ncol(Q_matrix_full) == 1 && nrow(Q_matrix_full) >= country_end_idx) { # Column vector
      X_country_native <- as.numeric(Q_matrix_full[country_start_idx:country_end_idx, 1])
    } else {
      warning(paste("Q_matrix for year", eora_data$year, "has unexpected dimensions (", 
                    nrow(Q_matrix_full), "x", ncol(Q_matrix_full), ") for country", country_iso3c, 
                    ". Will attempt to estimate total output X from Z sums if needed."))
    }
  }
  
  # If X_country_native couldn't be determined from Q_matrix, estimate from Z sums
  # This is a fallback and might be less accurate as it ignores final demand components in X.
  if (is.null(X_country_native)) {
    X_country_native <- rowSums(Z_country_native) # Sum of intermediate deliveries by sector
    # This is an approximation. Ideally X = rowSums(Z) + rowSums(FinalDemandMatrix)
    warning(paste("Estimating total output X for", country_iso3c, "year", eora_data$year, 
                  "from Z sums (row sums of intermediate transactions). Quality may be affected."))
  }
  
  # Apply PPP Scaling Factor to monetary flow data (Z and X)
  # Note: If ppp_scaling_factor is 1.0, these operations have no effect.
  Z_country_native_scaled <- Z_country_native * ppp_scaling_factor
  X_country_native_scaled <- X_country_native * ppp_scaling_factor
  
  if (length(X_country_native_scaled) != n_sectors_for_this_country) {
    warning(paste0("Scaled X_native length (", length(X_country_native_scaled), 
                   ") != n_sectors_label (", n_sectors_for_this_country, ") for ", 
                   country_iso3c, " yr ", eora_data$year))
    # Handle specific cases like 'ROW' if it appears with inconsistent dimensions
    if (length(X_country_native_scaled) > n_sectors_for_this_country && 
        n_sectors_for_this_country == 1 && country_iso3c == "ROW") {
      warning("ROW has 1 sector in labels, but X_scaled has more. Summing X_scaled.")
      X_country_native_scaled <- sum(X_country_native_scaled, na.rm = TRUE)
    } else {
      return(list(A_matrix = NULL, status = paste("X_native_scaled length mismatch for", country_iso3c)))
    }
  }
  
  # Handle non-positive or very small total outputs (X_country_native_scaled)
  # These can cause issues (division by zero) when calculating A_matrix
  if (any(X_country_native_scaled <= 1e-9, na.rm = TRUE)) {
    X_country_native_scaled[is.na(X_country_native_scaled)] <- 0 # Convert NAs to 0 before patching
    # Find mean of positive outputs to use for patching, or a small default if all are non-positive
    positive_X_mean <- mean(X_country_native_scaled[X_country_native_scaled > 1e-9], na.rm = TRUE)
    if (is.na(positive_X_mean) || !is.finite(positive_X_mean) || positive_X_mean <= 1e-9) {
      positive_X_mean <- 1 # Default small positive value if no valid mean
    }
    # Patch non-positive/small values with a fraction of the positive mean
    X_country_native_scaled[X_country_native_scaled <= 1e-9] <- positive_X_mean * 0.001 # Use a very small positive number
    warning(paste("Non-positive/small total output values (X_scaled) patched for", country_iso3c, 
                  "year", eora_data$year, "to avoid division by zero. This might affect A_matrix quality."))
  }
  
  country_sectors_native_labels <- eora_data$sector_labels$sector[country_start_idx:country_end_idx]
  if (length(country_sectors_native_labels) != n_sectors_for_this_country) {
    return(list(A_matrix = NULL, status = paste("Sector label length mismatch post-extraction for", country_iso3c, "yr", eora_data$year)))
  }
  
  # Handle 'ROW' if it's a single 'TOTAL' sector that cannot be meaningfully aggregated
  if (n_sectors_for_this_country == 1 && country_iso3c == "ROW" && 
      (length(country_sectors_native_labels) == 1 && toupper(country_sectors_native_labels) == "TOTAL")) {
    warning(paste("Skipping aggregation for ROW (single 'TOTAL' sector) for yr", eora_data$year))
    return(list(A_matrix = NULL, status = "ROW has single TOTAL sector, cannot aggregate"))
  }
  
  # Pass the SCALED Z and X matrices to aggregation function
  return(aggregate_to_sectors(Z_country_native_scaled, X_country_native_scaled, country_sectors_native_labels))
}

aggregate_to_sectors <- function(Z_native_scaled, X_native_scaled, native_sectors) {
  n_agg <- length(AGG_SECTOR_NAMES)
  # Initialize aggregated Z and X matrices/vectors with zeros
  Z_agg_scaled <- matrix(0, n_agg, n_agg, dimnames = list(AGG_SECTOR_NAMES, AGG_SECTOR_NAMES))
  X_agg_scaled <- stats::setNames(numeric(n_agg), AGG_SECTOR_NAMES)
  
  # Map native EORA sectors to the 8 aggregate sectors
  sector_mapping_vec <- sapply(native_sectors, map_eora_sector_to_8agg)
  valid_mapping_indices <- which(!is.na(sector_mapping_vec)) # Indices of native sectors that map to an aggregate sector
  
  if (length(valid_mapping_indices) == 0) {
    return(list(A_matrix = NULL, status = "No valid EORA sector maps for aggregation"))
  }
  
  # Aggregate Z_native_scaled and X_native_scaled
  for (i in valid_mapping_indices) { # Iterate over rows (supplying native sector)
    agg_row_sector <- sector_mapping_vec[i]
    if (!is.na(X_native_scaled[i])) {
      X_agg_scaled[agg_row_sector] <- X_agg_scaled[agg_row_sector] + X_native_scaled[i]
    }
    for (j in valid_mapping_indices) { # Iterate over columns (receiving native sector)
      agg_col_sector <- sector_mapping_vec[j]
      if (!is.na(Z_native_scaled[i, j])) {
        Z_agg_scaled[agg_row_sector, agg_col_sector] <- Z_agg_scaled[agg_row_sector, agg_col_sector] + Z_native_scaled[i, j]
      }
    }
  }
  
  # Calculate A_agg (technical coefficient matrix) from aggregated & scaled Z and X
  # A_ij = Z_ij / X_j. If Z and X were scaled by the same factor, A should be invariant.
  A_agg <- Z_agg_scaled # Initialize A_agg with Z_agg_scaled values
  for (j_col_idx in 1:n_agg) { # For each column (receiving sector j)
    if (X_agg_scaled[j_col_idx] > 1e-9) { # Avoid division by zero or very small numbers
      A_agg[, j_col_idx] <- Z_agg_scaled[, j_col_idx] / X_agg_scaled[j_col_idx]
    } else {
      A_agg[, j_col_idx] <- 0 # If total output of sector j is effectively zero, inputs to it are zero
    }
  }
  
  A_agg[!is.finite(A_agg)] <- 0 # Replace Inf/NaN with 0 (e.g., from 0/0)
  
  # Optional: Check and normalize columns of A_agg if sums are > 1 (economically implausible)
  # col_sums_A <- colSums(A_agg)
  # problematic_cols <- which(col_sums_A > 0.99) # Using 0.99 as a threshold for potential issues
  # if (length(problematic_cols) > 0) {
  #   warning(paste("Columns in A_agg matrix sum to >0.99 for sectors:", 
  #                 paste(AGG_SECTOR_NAMES[problematic_cols], collapse=", "), 
  #                 ". This might indicate issues with data or aggregation. Scaling them down."))
  #   for (col_idx in problematic_cols) {
  #     if (col_sums_A[col_idx] > 0) { # Avoid division by zero if sum is already zero
  #       A_agg[, col_idx] <- A_agg[, col_idx] * (0.98 / col_sums_A[col_idx]) # Scale to sum to 0.98
  #     }
  #   }
  # }
  
  return(list(A_matrix = A_agg, status = "Success", X_vector_agg_scaled = X_agg_scaled, Z_matrix_agg_scaled = Z_agg_scaled))
}

compute_io_measures_revised <- function(A_matrix, manuf_indices, serv_indices) {
  n <- nrow(A_matrix)
  I_mat <- diag(n)
  
  # Attempt to compute Leontief inverse L_inv = (I - A)^-1
  # Add small epsilon to diagonal of (I-A) if it's singular, to aid inversion
  L_inv <- tryCatch(
    solve(I_mat - A_matrix),
    error = function(e) {
      warning(paste("Matrix (I-A) is singular or near-singular. Attempting pseudo-inverse or regularization for L_inv. Error:", e$message))
      # Try adding a small value to diagonal of (I-A) for stability
      solve(I_mat - A_matrix + diag(1e-7, n)) 
    }
  )
  L_inv[L_inv < 0] <- 0 # Negative values in L_inv are economically not meaningful, set to 0.
  
  manuf_idx_valid <- manuf_indices[manuf_indices >= 1 & manuf_indices <= n]
  serv_idx_valid <- serv_indices[serv_indices >= 1 & serv_indices <= n]
  
  # Check if A_matrix is all zeros or L_inv is all zeros, which might indicate issues.
  if(all(A_matrix == 0)) warning("A_matrix is all zeros. IO measures will be trivial.")
  if(all(L_inv == 0)) warning("L_inv (Leontief Inverse) is all zeros. Multiplier effects will be zero.")
  
  return(list(
    A_matrix = A_matrix, # Original A_matrix passed in
    L_matrix = L_inv,    # This is the Leontief Inverse (I-A)^-1
    upstream_centrality = mean(colSums(A_matrix), na.rm = TRUE), # Avg column sums of A (impact of a sector on suppliers)
    downstream_centrality = mean(rowSums(A_matrix), na.rm = TRUE), # Avg row sums of A (impact of a sector on customers)
    multiplier_effect = mean(colSums(L_inv), na.rm = TRUE), # Avg column sums of L_inv (total output multiplier)
    manufacturing_intensity = if (length(manuf_idx_valid) > 0 && nrow(A_matrix) >= max(manuf_idx_valid, na.rm=TRUE)) mean(A_matrix[manuf_idx_valid, manuf_idx_valid, drop = FALSE], na.rm = TRUE) else NA_real_,
    services_intensity = if (length(serv_idx_valid) > 0 && nrow(A_matrix) >= max(serv_idx_valid, na.rm=TRUE)) mean(A_matrix[serv_idx_valid, serv_idx_valid, drop = FALSE], na.rm = TRUE) else NA_real_,
    economic_complexity = mean(colSums(A_matrix > 0.01, na.rm = TRUE), na.rm = TRUE) # Avg number of significant inputs per sector
  ))
}

run_scalar_rdd_analysis <- function(data_for_rdd, config_opts) {
  scalar_outcomes <- c("upstream_centrality", "downstream_centrality", "multiplier_effect", "manufacturing_intensity", "services_intensity", "economic_complexity")
  results <- list()
  
  if (nrow(data_for_rdd %||% data.frame()) == 0) {
    cat("No data for scalar RDD.\n")
    return(results)
  }
  
  for (outcome_name in scalar_outcomes) {
    min_obs_check <- config_opts$min_obs_per_side * 2
    if (!outcome_name %in% names(data_for_rdd)) {
      cat("Skip RDD for", outcome_name, "- not in data.\n")
      results[[outcome_name]] <- list(outcome = outcome_name, error = "Not in data")
      next
    }
    if (sum(!is.na(data_for_rdd[[outcome_name]])) < min_obs_check) {
      cat("Skip RDD for", outcome_name, "- insufficient non-NA obs (", sum(!is.na(data_for_rdd[[outcome_name]])), ").\n")
      results[[outcome_name]] <- list(outcome = outcome_name, error = paste("Insufficient non-NA obs:", sum(!is.na(data_for_rdd[[outcome_name]]))))
      next
    }
    
    cat("\nRunning RDD for:", outcome_name, "\n")
    current_data <- data_for_rdd %>%
      dplyr::select(dplyr::all_of(c(outcome_name, "running_var_centered", "treatment"))) %>% # Added treatment for N_left/N_right check
      dplyr::filter(!is.na(.data[[outcome_name]]), !is.na(running_var_centered))
    
    if (nrow(current_data) < min_obs_check) {
      cat("  Insufficient obs for", outcome_name, "with RV post-filter (", nrow(current_data), ").\n")
      results[[outcome_name]] <- list(outcome = outcome_name, error = "Insufficient data post-filter for RDD")
      next
    }
    # Check for min_obs_per_side after filtering for this specific outcome
    if (sum(current_data$treatment == 0) < config_opts$min_obs_per_side || 
        sum(current_data$treatment == 1) < config_opts$min_obs_per_side) {
      cat("  Insufficient obs on one side of cutoff for", outcome_name, 
          "(Left:", sum(current_data$treatment == 0), ", Right:", sum(current_data$treatment == 1), "). Min per side:", config_opts$min_obs_per_side, "\n")
      results[[outcome_name]] <- list(outcome = outcome_name, error = "Insufficient obs on one side of cutoff")
      next
    }
    
    results[[outcome_name]] <- tryCatch(
      {
        rd_model <- rdrobust::rdrobust(y = current_data[[outcome_name]], x = current_data$running_var_centered, c = 0)
        list(
          outcome = outcome_name,
          coefficient = rd_model$coef[1], # Conventional RD estimate
          se = rd_model$se[3],           # Robust SE
          pvalue = rd_model$pv[3],         # Robust p-value
          bandwidth = rd_model$bws[1, 1], # Bandwidth used for estimate
          n_left = rd_model$N[1],        # N obs left of cutoff within bandwidth
          n_right = rd_model$N[2],       # N obs right of cutoff within bandwidth
          rd_result_full = rd_model
        )
      },
      error = function(e) {
        cat("Error in rdrobust for", outcome_name, ":", e$message, "\n")
        list(outcome = outcome_name, error = e$message)
      }
    )
  }
  return(results)
}

run_frechet_tests <- function(matrices_A, matrices_L, named_running_vars, config_opts) {
  results <- list()
  if (!"frechesTest" %in% ls(envir = .GlobalEnv) && !requireNamespace("frechesTest", quietly = TRUE)) {
    # Try to source if not loaded and not a package, assuming it's a script.
    # This part is heuristic. Best practice is to ensure frechesTest is available via library() or source() explicitly at script start.
    if (file.exists("frechesTest.R")) { # A common name for a local script
      tryCatch(source("frechesTest.R"), error = function(e) warning("Failed to source frechesTest.R"))
    }
    if (!"frechesTest" %in% ls(envir = .GlobalEnv)) { # Check again
      warning("`frechesTest` function not available. Ensure it's loaded (e.g., from a package or sourced script). Skipping Fréchet tests.")
      return(results)
    }
  }
  
  if (is.null(names(named_running_vars)) || length(names(named_running_vars)) != length(named_running_vars)) {
    warning("`named_running_vars` for Fréchet test is not a named vector or names are missing. Skipping.")
    return(results)
  }
  
  # Using matrices_L which should be Leontief Inverse (I-A)^-1 from compute_io_measures_revised
  cat("\nFréchet Test: A-matrix (covariance space)...\n")
  results$A_covariance <- run_single_frechet_test(matrices_A, named_running_vars, "covariance", config_opts)
  
  cat("\nFréchet Test: L_inv-matrix (Leontief Inverse, covariance space)...\n")
  results$L_inv_covariance <- run_single_frechet_test(matrices_L, named_running_vars, "covariance", config_opts) # matrices_L is L_inv
  
  # For network space, we typically use Adjacency or Laplacians derived from Adjacency
  # The current convert_io_matrices_to_laplacians uses A as input.
  # If "network" space implies graph Laplacians, it's fine.
  cat("\nFréchet Test: A-matrix (as Adjacency for network space, then to Laplacian)...\n")
  results$A_network_laplacian <- run_single_frechet_test(matrices_A, named_running_vars, "network", config_opts)
  
  cat("\nFréchet Test: L_inv-matrix (as Adjacency for network space, then to Laplacian)...\n")
  results$L_inv_network_laplacian <- run_single_frechet_test(matrices_L, named_running_vars, "network", config_opts) # matrices_L is L_inv
  
  return(results)
}


run_single_frechet_test <- function(matrices_list_input, running_vars_input, metric_space_for_frechestest, config_opts) {
  # Filter out matrices or RVs with NA names (obs_ids)
  valid_matrix_names <- names(matrices_list_input)[!is.na(names(matrices_list_input))]
  valid_rv_names <- names(running_vars_input)[!is.na(names(running_vars_input))]
  
  matrices_list_input_clean <- matrices_list_input[valid_matrix_names]
  running_vars_input_clean <- running_vars_input[valid_rv_names]
  
  common_obs_ids <- intersect(names(matrices_list_input_clean), names(running_vars_input_clean))
  
  min_obs_frechet <- max(10, config_opts$cv_k_folds + 5, config_opts$min_obs_per_side * 2) # Ensure enough obs
  if (length(common_obs_ids) < min_obs_frechet) {
    error_msg <- paste("Insufficient common obs (", length(common_obs_ids), ") for Fréchet test (metric:", metric_space_for_frechestest, "). Min ", min_obs_frechet, " required.")
    cat(error_msg, "\n")
    return(list(p_value = NA, test_statistic = NA, error = error_msg, frechet_mean_C = NULL, frechet_mean_T = NULL, frechet_mean_diff = NULL))
  }
  
  # Align matrices and running variable by common observation IDs
  Y_obj_for_test_orig <- matrices_list_input_clean[common_obs_ids]
  X_scalar_for_test_orig <- running_vars_input_clean[common_obs_ids]
  
  # Filter out any NA/Inf in running variable or NULL/NA matrices before test
  valid_indices <- sapply(seq_along(Y_obj_for_test_orig), function(i) {
    !is.null(Y_obj_for_test_orig[[i]]) && 
      is.matrix(Y_obj_for_test_orig[[i]]) &&
      !any(is.na(Y_obj_for_test_orig[[i]])) && 
      !any(is.infinite(Y_obj_for_test_orig[[i]])) &&
      !is.na(X_scalar_for_test_orig[i]) &&
      is.finite(X_scalar_for_test_orig[i])
  })
  
  if (sum(valid_indices) < min_obs_frechet) {
    error_msg <- paste("Insufficient valid obs (", sum(valid_indices), ") after cleaning NA/Inf in Y or X for Fréchet test (metric:", metric_space_for_frechestest, "). Min ", min_obs_frechet, " required.")
    cat(error_msg, "\n")
    return(list(p_value = NA, test_statistic = NA, error = error_msg, frechet_mean_C = NULL, frechet_mean_T = NULL, frechet_mean_diff = NULL))
  }
  
  Y_obj_for_test <- Y_obj_for_test_orig[valid_indices]
  X_scalar_for_test <- X_scalar_for_test_orig[valid_indices]
  
  # If metric_space is "network", convert matrices to Laplacians
  # This conversion happens *after* initial filtering for NA/Inf in original matrices
  if (metric_space_for_frechestest == "network") {
    cat("  Converting matrices to Laplacians for network space test...\n")
    # Using symmetric normalized Laplacian as a common choice for network comparison
    Y_obj_laplacians <- convert_io_matrices_to_laplacians(Y_obj_for_test, method = "weighted", laplacian_type = "symmetric") 
    
    if (length(Y_obj_laplacians) < min_obs_frechet) {
      error_msg <- paste("Too few matrices (", length(Y_obj_laplacians), ") post-Laplacian conversion. Min ", min_obs_frechet, " required.")
      cat(error_msg, "\n")
      return(list(p_value = NA, test_statistic = NA, error = error_msg, frechet_mean_C = NULL, frechet_mean_T = NULL, frechet_mean_diff = NULL))
    }
    # Align X_scalar with the successfully converted Laplacians
    valid_laplacian_ids <- names(Y_obj_laplacians)
    X_scalar_for_test <- X_scalar_for_test[valid_laplacian_ids] # Assumes Y_obj_laplacians retains original names
    Y_obj_for_test <- Y_obj_laplacians # Use Laplacians for the test
  }
  
  # Final check on lengths
  if (length(Y_obj_for_test) != length(X_scalar_for_test)) {
    error_msg <- paste("CRITICAL: Y/X length mismatch before Fréchet test. Y:", length(Y_obj_for_test), "X:", length(X_scalar_for_test), "Metric:", metric_space_for_frechestest)
    cat(error_msg, "\n")
    return(list(p_value = NA, test_statistic = NA, error = error_msg, frechet_mean_C = NULL, frechet_mean_T = NULL, frechet_mean_diff = NULL))
  }
  
  # Check for sufficient observations on both sides of the cutoff
  num_left <- sum(X_scalar_for_test < 0, na.rm = TRUE)
  num_right <- sum(X_scalar_for_test >= 0, na.rm = TRUE)
  if (num_left < config_opts$min_obs_per_side || num_right < config_opts$min_obs_per_side) {
    error_msg <- paste0("Insufficient observations on one side of the cutoff for Fréchet test. ",
                        "Left: ", num_left, ", Right: ", num_right, 
                        ". Min per side required: ", config_opts$min_obs_per_side)
    cat(error_msg, "\n")
    return(list(p_value = NA, test_statistic = NA, error = error_msg, 
                frechet_mean_C = NULL, frechet_mean_T = NULL, frechet_mean_diff = NULL,
                n_left = num_left, n_right = num_right))
  }
  
  frechet_options_list <- list(metric = "frobenius") # Default metric for covariance space
  if (metric_space_for_frechestest == "network") { 
    # Options specific to network space in frechesTest package
    # These might need adjustment based on the specifics of the frechesTest implementation
    frechet_options_list$W_laplacian_bound <- 1.0 # Example option
    frechet_options_list$network_directed <- FALSE # Assuming undirected Laplacians
    # The metric for network space might be handled internally by frechesTest based on metric_space_type
  }
  
  cat("  Running frechesTest with", length(Y_obj_for_test), "valid obs for", metric_space_for_frechestest, "space.\n")
  
  # Determine function call based on availability (package vs. global)
  frechesTest_fn <- if (requireNamespace("frechesTest", quietly = TRUE)) {
    frechesTest::frechesTest
  } else if (exists("frechesTest", envir = .GlobalEnv)) {
    get("frechesTest", envir = .GlobalEnv)
  } else { # Should have been caught earlier
    stop("frechesTest function is not available.")
  }
  
  frechet_result_obj <- tryCatch(
    {
      frechesTest_fn( # Use determined function
        Y_obj = Y_obj_for_test, X_scalar = X_scalar_for_test, c_val = 0,
        metric_space_type = metric_space_for_frechestest, # "covariance" or "network"
        h_frechet = "CV", # Use Cross-Validation for bandwidth
        kernel_frechet_char = "epan", # Epanechnikov kernel
        frechet_optns = frechet_options_list,
        cv_K_folds = config_opts$cv_k_folds,
        cv_n_bw_candidates = config_opts$cv_n_bw_candidates, 
        verbose = FALSE # Set to TRUE for more detailed output from frechesTest
      )
    },
    error = function(e) {
      cat(" Error in frechesTest (metric:", metric_space_for_frechestest, "):", e$message, "\n")
      list(p_value = NA_real_, test_statistic = NA_real_, error = e$message, 
           frechet_mean_C = NULL, frechet_mean_T = NULL, frechet_mean_diff = NULL,
           n_left = num_left, n_right = num_right) # Include Ns even on error
    }
  )
  
  # Standardize output list
  return_list <- list(
    p_value = frechet_result_obj$p_value %||% frechet_result_obj$p.value %||% NA_real_, # frechesTest might use p.value
    test_statistic = frechet_result_obj$Tn %||% frechet_result_obj$test_statistic %||% NA_real_,
    bandwidth = frechet_result_obj$h_variance_used %||% frechet_result_obj$h_chosen %||% NA_real_,
    n_left = frechet_result_obj$N_left %||% frechet_result_obj$n_left %||% num_left %||% NA_integer_, # Use calculated Ns if available
    n_right = frechet_result_obj$N_right %||% frechet_result_obj$n_right %||% num_right %||% NA_integer_,
    error = frechet_result_obj$error %||% NA_character_,
    frechet_mean_C = NULL, 
    frechet_mean_T = NULL,
    frechet_mean_diff = NULL
  )
  
  # Assign Fréchet means if they exist in the result object
  # frechesTest might name these differently, e.g., l_hat_minus/plus or Y0_adj_mean_c/Y1_adj_mean_c
  if (!is.null(frechet_result_obj)) {
    if (!is.null(frechet_result_obj$l_hat_minus)) return_list$frechet_mean_C <- frechet_result_obj$l_hat_minus
    if (!is.null(frechet_result_obj$l_hat_plus)) return_list$frechet_mean_T <- frechet_result_obj$l_hat_plus
    
    # Fallback for other potential names from different versions/implementations
    if (is.null(return_list$frechet_mean_C) && !is.null(frechet_result_obj$Y0_adj_mean_c)) return_list$frechet_mean_C <- frechet_result_obj$Y0_adj_mean_c
    if (is.null(return_list$frechet_mean_T) && !is.null(frechet_result_obj$Y1_adj_mean_c)) return_list$frechet_mean_T <- frechet_result_obj$Y1_adj_mean_c
    
    if (is.null(return_list$frechet_mean_C) && !is.null(frechet_result_obj$frechet_mean_C)) return_list$frechet_mean_C <- frechet_result_obj$frechet_mean_C # if already named so
    if (is.null(return_list$frechet_mean_T) && !is.null(frechet_result_obj$frechet_mean_T)) return_list$frechet_mean_T <- frechet_result_obj$frechet_mean_T
    
    # Calculate difference: Treated - Control
    if (!is.null(return_list$frechet_mean_T) && !is.null(return_list$frechet_mean_C)) {
      if (is.matrix(return_list$frechet_mean_T) && is.matrix(return_list$frechet_mean_C) && 
          all(dim(return_list$frechet_mean_T) == dim(return_list$frechet_mean_C))) {
        return_list$frechet_mean_diff <- return_list$frechet_mean_T - return_list$frechet_mean_C
      } else {
        warning("Fréchet means (T,C) are not conformable matrices or not matrices for test: ", metric_space_for_frechestest, ". Cannot compute difference.")
      }
    }
  }
  
  # ————————————————————————————————————————————
  # if we're in "network" mode, pull the Laplacian means back to Leontief-inverse means
  if (metric_space_for_frechestest == "network" &&
      !is.null(return_list$frechet_mean_C) &&
      !is.null(return_list$frechet_mean_T)) {
    
    Lbar_C <- return_list$frechet_mean_C
    Lbar_T <- return_list$frechet_mean_T
    
    # reconstruct the adjacency (i.e. the Leontief inverse)
    Wbar_C <- diag(diag(Lbar_C)) - Lbar_C
    Wbar_T <- diag(diag(Lbar_T)) - Lbar_T
    
    # overwrite the Frechet‐mean slots with the adjacency means
    return_list$frechet_mean_C    <- Wbar_C
    return_list$frechet_mean_T    <- Wbar_T
    return_list$frechet_mean_diff <- Wbar_T - Wbar_C
  }
  # ————————————————————————————————————————————
  
  return(return_list)
}

run_analysis_with_io_lags <- function(wb_data, config_opts, base_rv_year, io_outcome_lags = c(1, 2, 3)) {
  cat("\n=== Running Analysis for RV Yr", base_rv_year, "with I-O Lags:", paste(io_outcome_lags, collapse = ", "), "===\n")
  results_by_io_lag <- list()
  
  temp_config_for_rv_prep <- config_opts 
  temp_config_for_rv_prep$analysis_mode <- "single_year"
  temp_config_for_rv_prep$single_year_choice <- base_rv_year
  
  base_data_for_rv_setup <- prepare_analysis_data(wb_data$countries_data, wb_data$thresholds, temp_config_for_rv_prep)
  
  if (is.null(base_data_for_rv_setup) || nrow(base_data_for_rv_setup %||% data.frame()) < (config_opts$min_obs_per_side * 2)) {
    cat("  Insufficient base data for RV year", base_rv_year, ". Aborting lag analysis.\n")
    return(results_by_io_lag)
  }
  cat("  Base data for RV year", base_rv_year, "prepared initially:", nrow(base_data_for_rv_setup), "obs.\n")
  
  for (io_lag in io_outcome_lags) {
    current_io_outcome_year <- base_rv_year + io_lag 
    cat("\n  Processing I-O Lag +", io_lag, "(I-O Year:", current_io_outcome_year, ")\n")
    
    if (current_io_outcome_year > config_opts$max_io_year_available) {
      cat("    Skipping I-O year", current_io_outcome_year, "> max available", config_opts$max_io_year_available, "\n")
      next
    }
    
    analysis_data_for_this_io_lag_step1 <- base_data_for_rv_setup %>%
      dplyr::mutate(
        io_target_year = current_io_outcome_year, 
        obs_id = create_obs_id(iso3c, io_target_year) 
      )
    
    # --- DEBUG POINT FOR LAGS: Before selecting out old PPP columns ---
    cat("\n--- DEBUG LAGS (run_analysis_with_io_lags) - STEP 1: After initial mutate for io_target_year ---\n")
    cat("Names in analysis_data_for_this_io_lag_step1:\n")
    print(names(analysis_data_for_this_io_lag_step1))
    cat("--- END DEBUG LAGS STEP 1 ---\n\n")
    
    analysis_data_for_this_io_lag_step2 <- analysis_data_for_this_io_lag_step1 %>%
      dplyr::select(-any_of(c("ppp_scaling_factor", "ppp_lcu_per_intl_dollar", "exch_lcu_per_usd", 
                              "gdp_ppp_curr_intl", "gdp_ppp_const_2017_intl")))
    
    # --- DEBUG POINT FOR LAGS: After selecting out old PPP columns ---
    cat("\n--- DEBUG LAGS (run_analysis_with_io_lags) - STEP 2: After selecting out old PPP columns ---\n")
    cat("Names in analysis_data_for_this_io_lag_step2 (should NOT have ppp_lcu_per_intl_dollar etc.):\n")
    print(names(analysis_data_for_this_io_lag_step2))
    cat("--- END DEBUG LAGS STEP 2 ---\n\n")
    
    # --- DEBUG POINT FOR LAGS: Check wb_data$countries_data just before it's used in join ---
    cat("\n--- DEBUG LAGS (run_analysis_with_io_lags) - PRE-JOIN: Names in wb_data$countries_data ---\n")
    print(names(wb_data$countries_data))
    cat("Structure of wb_data$countries_data (relevant columns):\n")
    if(all(c("iso3c", "year", "ppp_lcu_per_intl_dollar", "exch_lcu_per_usd") %in% names(wb_data$countries_data))){
      print(str(wb_data$countries_data %>% dplyr::select(iso3c, year, ppp_lcu_per_intl_dollar, exch_lcu_per_usd), list.len=5))
    } else {
      cat("CRITICAL FOR LAGS: One of the required PPP columns is missing from wb_data$countries_data!\n")
    }
    cat("--- END DEBUG LAGS PRE-JOIN ---\n\n")
    
    # This is the right-hand side of the join. If this select fails, it's because
    # wb_data$countries_data doesn't have one of these columns.
    ppp_data_for_join_rhs <- wb_data$countries_data %>% 
      dplyr::select(iso3c, year, ppp_lcu_per_intl_dollar, exch_lcu_per_usd, 
                    gdp_ppp_curr_intl, gdp_ppp_const_2017_intl) %>% 
      dplyr::rename(data_year_for_ppp = year)
    
    # --- DEBUG POINT FOR LAGS: After creating RHS for join ---
    cat("\n--- DEBUG LAGS (run_analysis_with_io_lags) - RHS FOR JOIN: Names in ppp_data_for_join_rhs ---\n")
    print(names(ppp_data_for_join_rhs))
    cat("--- END DEBUG LAGS RHS FOR JOIN ---\n\n")
    
    analysis_data_for_this_io_lag_step3 <- analysis_data_for_this_io_lag_step2 %>%
      dplyr::left_join(
        ppp_data_for_join_rhs, # Use the pre-created RHS
        by = c("iso3c", "io_target_year" = "data_year_for_ppp")
      )
    
    # --- DEBUG POINT FOR LAGS: After the join ---
    cat("\n--- DEBUG LAGS (run_analysis_with_io_lags) - STEP 3: After left_join ---\n")
    cat("Names in analysis_data_for_this_io_lag_step3 (should now have fresh PPP columns, no suffixes needed if step 2 worked):\n")
    print(names(analysis_data_for_this_io_lag_step3))
    if (!all(c("ppp_lcu_per_intl_dollar", "exch_lcu_per_usd") %in% names(analysis_data_for_this_io_lag_step3))) {
      cat("CRITICAL FOR LAGS JOIN: Fresh PPP columns NOT found after join! Check join conditions and RHS data.\n")
    }
    cat("--- END DEBUG LAGS STEP 3 ---\n\n")
    
    analysis_data_for_this_io_lag <- analysis_data_for_this_io_lag_step3 %>%
      dplyr::mutate(
        temp_pli_usd_per_intl_dollar = as.numeric(ppp_lcu_per_intl_dollar) / as.numeric(exch_lcu_per_usd),
        temp_deflator_intl_dollar_base2017 = as.numeric(gdp_ppp_curr_intl) / as.numeric(gdp_ppp_const_2017_intl),
        initial_ppp_scaling_factor = (1 / temp_pli_usd_per_intl_dollar) / temp_deflator_intl_dollar_base2017,
        ppp_scaling_factor = ifelse(is.na(initial_ppp_scaling_factor) | !is.finite(initial_ppp_scaling_factor) | initial_ppp_scaling_factor <= 0,
                                    1.0, initial_ppp_scaling_factor)
      ) %>%
      dplyr::select(-dplyr::starts_with("temp_"), -initial_ppp_scaling_factor) %>%
      dplyr::distinct(obs_id, .keep_all = TRUE)
    
    num_defaulted_lag <- sum(analysis_data_for_this_io_lag$ppp_scaling_factor == 1.0 & 
                               (is.na(analysis_data_for_this_io_lag$ppp_lcu_per_intl_dollar) | 
                                  !is.finite(analysis_data_for_this_io_lag$ppp_lcu_per_intl_dollar / analysis_data_for_this_io_lag$exch_lcu_per_usd) |
                                  (analysis_data_for_this_io_lag$ppp_lcu_per_intl_dollar / analysis_data_for_this_io_lag$exch_lcu_per_usd) <=0 ), na.rm=TRUE) 
    cat("    For I-O lag +", io_lag, "(IO year ", current_io_outcome_year, "), ", 
        nrow(analysis_data_for_this_io_lag), " obs. ",
        num_defaulted_lag, " defaulted to PPP factor 1.0.\n")
    
    if (nrow(analysis_data_for_this_io_lag) < (config_opts$min_obs_per_side * 2)) {
      cat("    Insufficient data after setting I-O year to", current_io_outcome_year, 
          " (", nrow(analysis_data_for_this_io_lag), " obs). Skipping this lag.\n")
      next
    }
    
    io_data_for_this_lag <- process_eora_data(analysis_data_for_this_io_lag, config_opts)
    
    if (is.null(io_data_for_this_lag) || nrow(io_data_for_this_lag$data %||% data.frame()) < config_opts$min_obs_per_side * 2) {
      cat("    Insufficient I-O data after processing for I-O year", current_io_outcome_year, ". Skipping this lag.\n")
      next
    }
    
    named_rv_for_this_lag <- stats::setNames(io_data_for_this_lag$data$running_var_centered, io_data_for_this_lag$data$obs_id)
    results_by_io_lag[[paste0("rv_", base_rv_year, "_io_lag_", io_lag)]] <- list(
      rv_year = base_rv_year,
      io_lag_applied = io_lag,
      io_year_used = current_io_outcome_year,
      n_obs_after_io_processing = nrow(io_data_for_this_lag$data),
      rdd_results = run_scalar_rdd_analysis(io_data_for_this_lag$data, config_opts),
      frechet_results = run_frechet_tests(
        io_data_for_this_lag$matrices_A,
        io_data_for_this_lag$matrices_L, 
        named_rv_for_this_lag, config_opts
      )
    )
  }
  return(results_by_io_lag)
}

generate_all_outputs <- function(rdd_results, frechet_results, lagged_io_results, final_analysis_data, config_opts) {
  output_dir_path <- config_opts$output_dir
  raw_results_dir <- file.path(output_dir_path, "raw_results")
  plots_frechet_dir <- file.path(output_dir_path, "plots", "frechet_diff") # Defined earlier in setup_directories
  # if (!dir.exists(raw_results_dir)) dir.create(raw_results_dir, recursive = TRUE) # Done by setup_directories
  # if (!dir.exists(plots_frechet_dir)) dir.create(plots_frechet_dir, recursive = TRUE) # Done by setup_directories
  
  cat("\nCreating summary tables...\n")
  scalar_summary_df <- create_scalar_summary(rdd_results)
  utils::write.csv(scalar_summary_df, file.path(output_dir_path, "tables", "scalar_rdd_summary.csv"), row.names = FALSE)
  frechet_summary_df <- create_frechet_summary(frechet_results)
  utils::write.csv(frechet_summary_df, file.path(output_dir_path, "tables", "frechet_test_summary.csv"), row.names = FALSE)
  lagged_summary_df <- NULL
  if (length(lagged_io_results) > 0) {
    lagged_summary_df <- create_lagged_summary(lagged_io_results)
    utils::write.csv(lagged_summary_df, file.path(output_dir_path, "tables", "lagged_analysis_summary.csv"), row.names = FALSE)
  }
  
  cat("\nSaving full Fréchet mean estimates (main analysis)...\n")
  for (test_name in names(frechet_results)) {
    result_item <- frechet_results[[test_name]]
    if (!is.null(result_item$frechet_mean_C) && !is.null(result_item$frechet_mean_T)) {
      frechet_means_to_save <- list(mean_control = result_item$frechet_mean_C, 
                                    mean_treated = result_item$frechet_mean_T, 
                                    difference = result_item$frechet_mean_diff) # Should be T-C
      saveRDS(frechet_means_to_save, file = file.path(raw_results_dir, paste0("frechet_means_main_", test_name, ".Rds")))
      cat("  Saved Fréchet means for main test:", test_name, "\n")
    } else if (is.null(result_item$error) || is.na(result_item$error)) {
      warning("Fréchet means not found or not calculable in result object for main test: ", test_name)
    }
  }
  
  cat("\nSaving full Fréchet mean estimates (lagged analysis)...\n")
  if (length(lagged_io_results) > 0) {
    for (lag_result_key in names(lagged_io_results)) {
      lag_item <- lagged_io_results[[lag_result_key]]
      if (!is.null(lag_item$frechet_results)) {
        for (test_name in names(lag_item$frechet_results)) {
          result_item <- lag_item$frechet_results[[test_name]]
          if (!is.null(result_item$frechet_mean_C) && !is.null(result_item$frechet_mean_T)) {
            frechet_means_to_save <- list(mean_control = result_item$frechet_mean_C, 
                                          mean_treated = result_item$frechet_mean_T, 
                                          difference = result_item$frechet_mean_diff) # Should be T-C
            saveRDS(frechet_means_to_save, file = file.path(raw_results_dir, paste0("frechet_means_lagged_", lag_result_key, "_", test_name, ".Rds")))
            cat("  Saved Fréchet means for lagged test:", lag_result_key, "_", test_name, "\n")
          } else if (is.null(result_item$error) || is.na(result_item$error)) {
            warning("Fréchet means not found or not calculable in result object for lagged test: ", lag_result_key, "_", test_name)
          }
        }
      }
    }
  }
  
  cat("\nCreating plots...\n")
  if (nrow(final_analysis_data %||% data.frame()) > 0) {
    scalar_outcomes_to_plot <- names(rdd_results)
    for (outcome_plot_name in scalar_outcomes_to_plot) {
      # Check if result exists and is not an error object (i.e., has a coefficient)
      if (!is.null(rdd_results[[outcome_plot_name]]) && !is.null(rdd_results[[outcome_plot_name]]$coefficient) && !is.na(rdd_results[[outcome_plot_name]]$coefficient)) {
        rdd_plot_obj <- plot_rdd(final_analysis_data, outcome_plot_name, 
                                 bandwidth = rdd_results[[outcome_plot_name]]$bandwidth %||% NULL, 
                                 config_opts = config_opts)
        ggplot2::ggsave(
          filename = file.path(output_dir_path, "plots", paste0("rdd_plot_", outcome_plot_name, ".png")),
          plot = rdd_plot_obj, width = 10, height = 8, dpi = 300
        )
      } else {
        cat("Skipping plot for", outcome_plot_name, "due to missing, NA, or error in RDD result.\n")
      }
    }
    
    manipulation_test_result <- check_manipulation(final_analysis_data, config_opts)
    if (!is.null(manipulation_test_result) && inherits(manipulation_test_result, "rddensity")) {
      grDevices::png(file.path(output_dir_path, "plots", "manipulation_check_density.png"), width = 800, height = 600)
      # Need to handle if plot.rddensity is not found or if it's an S3 method.
      # Assuming rddensity package is loaded and plot method is available.
      tryCatch(plot(manipulation_test_result), # Generic plot might dispatch to rddensity::plot.rddensity
               error = function(e) {
                 tryCatch(rddensity::plot.rddensity(manipulation_test_result), # Explicit call
                          error = function(e2) cat("Failed to plot rddensity result:", e2$message))
               })
      graphics::title(main = "Density of Running Variable (Manipulation Check)")
      grDevices::dev.off()
      cat("Manipulation test (rddensity) p-value (H0: no manipulation):", manipulation_test_result$p_jk %||% NA_real_, "\n")
    } else if (!is.null(manipulation_test_result)) { # Result exists but not 'rddensity' class
      warning("Manipulation test result was not NULL but not of class 'rddensity'. Skipping plot.")
    }
  } else {
    cat("Skipping plot generation: `final_analysis_data` is empty or NULL.\n")
  }
  
  cat("\nCreating Fréchet mean difference plots (main analysis)...\n")
  for (test_name in names(frechet_results)) {
    result_item <- frechet_results[[test_name]]
    if (!is.null(result_item$frechet_mean_C) && !is.null(result_item$frechet_mean_T)) {
      mean_C_plot <- result_item$frechet_mean_C
      mean_T_plot <- result_item$frechet_mean_T
      
      if (is.matrix(mean_C_plot) && is.matrix(mean_T_plot) &&
          nrow(mean_C_plot) == length(AGG_SECTOR_NAMES) && ncol(mean_C_plot) == length(AGG_SECTOR_NAMES) &&
          nrow(mean_T_plot) == length(AGG_SECTOR_NAMES) && ncol(mean_T_plot) == length(AGG_SECTOR_NAMES)) {
        
        # at top of your script, define:
        agg_labels_short <- c(
          AGR     = "Agri. & Fish.",
          MIN     = "Mining",
          MAN     = "Manufacturing",
          UTL     = "Utilities",
          CON     = "Construction",
          TRD     = "Trade & Hosp.",
          FIN_PRO = "Finance & Business",
          PUB_OTH = "Public & Other"
        )
        
        # then inside the plotting loop, before plot_frechet_mean_difference():
        short_names <- unname(agg_labels_short[AGG_SECTOR_NAMES])
        dimnames(mean_C_plot) <- list(short_names, short_names)
        dimnames(mean_T_plot) <- list(short_names, short_names)
        
        current_plot_title_suffix <- paste0("(", stringr::str_to_title(gsub("_", " ", test_name)), ", Main Analysis)")
        
        frechet_diff_plot <- plot_frechet_mean_difference(
          left_mean = mean_C_plot, # Control
          right_mean = mean_T_plot, # Treated
          title_suffix = current_plot_title_suffix
        )
        
        if (!is.null(frechet_diff_plot) && inherits(frechet_diff_plot, "ggplot")) {
          plot_filename <- paste0("frechet_diff_main_", test_name, ".png")
          ggplot2::ggsave(
            filename = file.path(plots_frechet_dir, plot_filename),
            plot = frechet_diff_plot,
            width = 8, height = 7, dpi = 300 
          )
          cat("  Saved Fréchet difference plot:", plot_filename, "\n")
        } else {
          warning(paste("Failed to generate Fréchet difference plot for main test:", test_name))
        }
      } else {
        warning(paste("Fréchet means for main test", test_name, "are not matrices, dimensions mismatch AGG_SECTOR_NAMES, or null. Skipping plot."))
      }
    }
  }
  
  cat("\nCreating Fréchet mean difference plots (lagged analysis)...\n")
  if (length(lagged_io_results) > 0) {
    for (lag_result_key in names(lagged_io_results)) {
      lag_item <- lagged_io_results[[lag_result_key]]
      if (!is.null(lag_item$frechet_results)) {
        for (test_name in names(lag_item$frechet_results)) {
          result_item <- lag_item$frechet_results[[test_name]]
          if (!is.null(result_item$frechet_mean_C) && !is.null(result_item$frechet_mean_T)) {
            mean_C_plot <- result_item$frechet_mean_C
            mean_T_plot <- result_item$frechet_mean_T
            
            if (is.matrix(mean_C_plot) && is.matrix(mean_T_plot) &&
                nrow(mean_C_plot) == length(AGG_SECTOR_NAMES) && ncol(mean_C_plot) == length(AGG_SECTOR_NAMES) &&
                nrow(mean_T_plot) == length(AGG_SECTOR_NAMES) && ncol(mean_T_plot) == length(AGG_SECTOR_NAMES)) {
              
              # at top of your script, define:
              agg_labels_short <- c(
                AGR     = "Agri. & Fish.",
                MIN     = "Mining",
                MAN     = "Manufacturing",
                UTL     = "Utilities",
                CON     = "Construction",
                TRD     = "Trade & Hosp.",
                FIN_PRO = "Finance & Business",
                PUB_OTH = "Public & Other"
              )
              
              # then inside the plotting loop, before plot_frechet_mean_difference():
              short_names <- unname(agg_labels_short[AGG_SECTOR_NAMES])
              dimnames(mean_C_plot) <- list(short_names, short_names)
              dimnames(mean_T_plot) <- list(short_names, short_names)
              
              current_plot_title_suffix <- paste0("(", stringr::str_to_title(gsub("_", " ", test_name)), ", ", lag_result_key, ")")
              
              frechet_diff_plot_lagged <- plot_frechet_mean_difference(
                left_mean = mean_C_plot, # Control
                right_mean = mean_T_plot, # Treated
                title_suffix = current_plot_title_suffix
              )
              
              if (!is.null(frechet_diff_plot_lagged) && inherits(frechet_diff_plot_lagged, "ggplot")) {
                plot_filename <- paste0("frechet_diff_lagged_", lag_result_key, "_", test_name, ".png")
                ggplot2::ggsave(
                  filename = file.path(plots_frechet_dir, plot_filename),
                  plot = frechet_diff_plot_lagged,
                  width = 8, height = 7, dpi = 300 
                )
                cat("  Saved Fréchet difference plot:", plot_filename, "\n")
              } else {
                warning(paste("Failed to generate Fréchet difference plot for lagged test:", lag_result_key, "_", test_name))
              }
            } else {
              warning(paste("Fréchet means for lagged test", lag_result_key, "_", test_name, "are not matrices, dimensions mismatch AGG_SECTOR_NAMES, or null. Skipping plot."))
            }
          }
        }
      }
    }
    
    # ———————————————————————————————————————————————————————————————————
    # Build and write a LaTeX table of the 2-year-lag scalar RDD estimates
    # ———————————————————————————————————————————————————————————————————
    # ———————————————————————————————————————————————————————————————————
    # Write a plain tabular (no table env) for the 2-year-lag scalar RDDs
    # ———————————————————————————————————————————————————————————————————
    if (exists("lagged_summary_df") && nrow(lagged_summary_df) > 0) {
      library(dplyr)
      
      # 1) pull out just the IO_Lag==2, scalar RDD rows
      df2 <- lagged_summary_df %>%
        filter(IO_Lag == 2, ResultType == "RDD Scalar") %>%
        select(Outcome, Estimate, Std.Error) %>%
        arrange(Outcome)
      
      # 2) get the universal N obs
      nobs <- unique(lagged_summary_df$`N Obs (Analysis)`[lagged_summary_df$IO_Lag == 2])
      if(length(nobs)>1) nobs <- nobs[1]
      
      # 3) build the lines of the tabular
      lines <- c(
        "\\begin{tabular}[t]{lr}",
        "\\toprule",
        "Outcome & Coefficient \\\\",
        "\\midrule"
      )
      
      for(i in seq_len(nrow(df2))) {
        o   <- df2$Outcome[i]
        est <- sprintf("%.3f", df2$Estimate[i])
        se  <- sprintf("(%.3f)", df2$Std.Error[i])
        
        lines <- c(lines,
                   paste0(o, " & ", est, " \\\\"),
                   paste0("  & ", se, " \\\\")
        )
      }
      
      lines <- c(lines,
                 "\\midrule",
                 paste0("\\multicolumn{1}{l}{Observations} & ",
                        "\\multicolumn{1}{r}{", nobs, "} \\\\"),
                 "\\bottomrule",
                 "\\end{tabular}"
      )
      
      writeLines(lines,
                 file.path(config_opts$output_dir, "tables", "lag2_scalar_rdd.tex"))
    }
  }
  
  # Analyze top drops and increases in Leontief inverse matrix
  cat("\n=== TOP CHANGES IN LEONTIEF INVERSE MATRIX (MAIN ANALYSIS) ===\n")
  leontief_changes_table <- analyze_leontief_changes(frechet_results)
  if (!is.null(leontief_changes_table) && nrow(leontief_changes_table) > 0) {
    print(leontief_changes_table)
    utils::write.csv(leontief_changes_table, 
                     file.path(output_dir_path, "tables", "leontief_top_changes_main.csv"), 
                     row.names = FALSE)
    cat("\nLeontief inverse matrix top changes saved to: leontief_top_changes_main.csv\n")
  }
  
  # Analyze top drops and increases for lagged analysis
  if (length(lagged_io_results) > 0) {
    cat("\n=== TOP CHANGES IN LEONTIEF INVERSE MATRIX (LAGGED ANALYSIS) ===\n")
    for (lag_result_key in names(lagged_io_results)) {
      lag_item <- lagged_io_results[[lag_result_key]]
      if (!is.null(lag_item$frechet_results)) {
        cat("\nAnalyzing changes for", lag_result_key, ":\n")
        leontief_changes_lagged <- analyze_leontief_changes(lag_item$frechet_results)
        if (!is.null(leontief_changes_lagged) && nrow(leontief_changes_lagged) > 0) {
          print(leontief_changes_lagged)
          utils::write.csv(leontief_changes_lagged, 
                           file.path(output_dir_path, "tables", paste0("leontief_top_changes_", lag_result_key, ".csv")), 
                           row.names = FALSE)
          cat("\nSaved to:", paste0("leontief_top_changes_", lag_result_key, ".csv"), "\n")
        }
      }
    }
  }
  
  cat("\n=== CONSOLE ANALYSIS SUMMARY ===\n")
  cat("Scalar RDD Results (Main Analysis):\n")
  print(scalar_summary_df)
  cat("\nFréchet Test Results (Main Analysis):\n")
  print(frechet_summary_df)
  if (!is.null(lagged_summary_df) && nrow(lagged_summary_df %||% data.frame()) > 0) {
    cat("\nLagged Analysis (Fixed RV, Lagged IO) Summary:\n")
    print(lagged_summary_df)
  }
}

analyze_leontief_changes <- function(frechet_results_list) {
  # Extract Leontief inverse matrix Fréchet means from results
  # Looking for L_inv_network_laplacian or L_inv_covariance results
  
  leontief_results <- frechet_results_list[grep("^L_inv", names(frechet_results_list))]
  
  if (length(leontief_results) == 0) {
    cat("No Leontief inverse matrix results found in Fréchet test results.\n")
    return(NULL)
  }
  
  # Find the first result with valid Fréchet means
  valid_result <- NULL
  test_name_used <- NULL
  
  for (test_name in names(leontief_results)) {
    result <- leontief_results[[test_name]]
    if (!is.null(result$frechet_mean_C) && !is.null(result$frechet_mean_T) &&
        is.matrix(result$frechet_mean_C) && is.matrix(result$frechet_mean_T)) {
      valid_result <- result
      test_name_used <- test_name
      break
    }
  }
  
  if (is.null(valid_result)) {
    cat("No valid Fréchet means found for Leontief inverse matrices.\n")
    return(NULL)
  }
  
  cat("Analyzing changes from test:", test_name_used, "\n")
  
  # Get control and treated means
  mean_C <- valid_result$frechet_mean_C
  mean_T <- valid_result$frechet_mean_T
  
  # Check dimensions match
  if (!all(dim(mean_C) == dim(mean_T))) {
    cat("Dimension mismatch between control and treated Fréchet means.\n")
    return(NULL)
  }
  
  # Calculate differences (Treated - Control)
  diff_matrix <- mean_T - mean_C
  
  # Check if matrices are symmetric (within numerical tolerance)
  is_symmetric <- function(mat, tol = 1e-10) {
    all(abs(mat - t(mat)) < tol)
  }
  
  symmetric_matrix <- is_symmetric(mean_C) && is_symmetric(mean_T)
  if (symmetric_matrix) {
    cat("Detected symmetric Leontief matrices. Using upper triangle only to avoid duplicates.\n")
  } else {
    cat("Matrices are not symmetric. Will process all elements.\n")
  }
  
  # Get sector names if available
  sector_names <- if (!is.null(rownames(mean_C))) {
    rownames(mean_C)
  } else if (nrow(mean_C) == length(AGG_SECTOR_NAMES)) {
    AGG_SECTOR_NAMES
  } else {
    paste0("Sector_", 1:nrow(mean_C))
  }
  
  # Create a data frame of all industry pairs with their changes
  all_changes <- data.frame()
  
  if (symmetric_matrix) {
    # For symmetric matrices, only consider upper triangle (excluding diagonal)
    for (i in 1:nrow(diff_matrix)) {
      for (j in i:ncol(diff_matrix)) {
        if (i != j) {  # Exclude diagonal elements
          all_changes <- rbind(all_changes, data.frame(
            Sector_1 = sector_names[i],
            Sector_2 = sector_names[j],
            Change = diff_matrix[i, j],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  } else {
    # For non-symmetric matrices, consider all off-diagonal elements
    for (i in 1:nrow(diff_matrix)) {
      for (j in 1:ncol(diff_matrix)) {
        if (i != j) {  # Exclude diagonal elements
          all_changes <- rbind(all_changes, data.frame(
            Sector_1 = sector_names[i],
            Sector_2 = sector_names[j],
            Change = diff_matrix[i, j],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Sort by absolute change magnitude
  all_changes <- all_changes[order(abs(all_changes$Change), decreasing = TRUE), ]
  
  # Get top 5 increases (positive changes)
  top_increases <- head(all_changes[all_changes$Change > 0, ], 5)
  top_increases$Change_Type <- "Increase"
  
  # Get top 5 decreases (negative changes)
  top_decreases <- head(all_changes[all_changes$Change < 0, ], 5)
  top_decreases$Change_Type <- "Decrease"
  
  # Combine results
  top_changes <- rbind(top_increases, top_decreases)
  
  # Sort by change type and magnitude
  top_changes <- top_changes[order(top_changes$Change_Type, -abs(top_changes$Change)), ]
  
  # Format the output table
  top_changes$Rank <- ave(abs(top_changes$Change), 
                          top_changes$Change_Type, 
                          FUN = function(x) rank(-x, ties.method = "first"))
  
  # Reorder columns and format
  top_changes <- top_changes[, c("Rank", "Change_Type", "Sector_1", "Sector_2", "Change")]
  top_changes$Change <- round(top_changes$Change, 6)
  
  # Add interpretation column
  top_changes$Interpretation <- paste0(
    "Linkage between ", top_changes$Sector_1, 
    " and ", top_changes$Sector_2,
    ifelse(top_changes$Change > 0, " increased by ", " decreased by "),
    format(abs(top_changes$Change), scientific = FALSE)
  )
  
  return(top_changes)
}

create_scalar_summary <- function(rdd_results_list) {
  if (length(rdd_results_list) == 0) {
    return(data.frame()) # Return empty data.frame if list is empty
  }
  dplyr::bind_rows(lapply(rdd_results_list, function(x) {
    # Ensure x is a list and has expected components, otherwise fill with NA
    if (is.list(x)) {
      data.frame(
        Outcome = stringr::str_to_title(stringr::str_replace_all(x$outcome %||% "N/A", "_", " ")),
        Coefficient = round(x$coefficient %||% NA_real_, 4),
        `Std Error` = round(x$se %||% NA_real_, 4),
        `P-value` = round(x$pvalue %||% NA_real_, 4),
        Bandwidth = round(x$bandwidth %||% NA_real_, 2),
        `N Left` = x$n_left %||% NA_integer_,
        `N Right` = x$n_right %||% NA_integer_,
        Error = shorten_strings(x$error %||% NA_character_, 40), # Keep error short
        check.names = FALSE
      )
    } else {
      # If x is not a list (e.g. an error message string itself, though unlikely with current structure)
      data.frame(Outcome = "Error in result format", Coefficient=NA, `Std Error`=NA, `P-value`=NA, Bandwidth=NA, `N Left`=NA, `N Right`=NA, Error=shorten_strings(as.character(x),40), check.names=FALSE)
    }
  }))
}

create_frechet_summary <- function(frechet_results_list) {
  if (length(frechet_results_list) == 0) {
    return(data.frame())
  }
  dplyr::bind_rows(lapply(names(frechet_results_list), function(test_name_key) {
    result_item <- frechet_results_list[[test_name_key]]
    # Default matrix type and metric space if parsing fails
    matrix_type_parsed <- "Unknown Matrix"
    metric_space_parsed <- "Unknown Space"
    
    parts <- stringr::str_split(test_name_key, "_")[[1]]
    if (length(parts) >= 2) {
      matrix_type_parsed <- paste0(parts[1], ifelse(length(parts)>2 && parts[2]!="covariance" && parts[2]!="network", paste0("_",parts[2]), ""), " matrix") # e.g. L_inv matrix
      metric_space_parsed <- parts[length(parts)] # Last part is metric space
      if (metric_space_parsed=="laplacian") metric_space_parsed <- "network (Laplacian)"
    }
    
    if (is.list(result_item)) {
      data.frame(
        `Test Name` = test_name_key, # Keep original test name
        `Matrix Type` = matrix_type_parsed,
        `Metric Space` = metric_space_parsed,
        `Test Statistic` = round(result_item$test_statistic %||% NA_real_, 4),
        `P-value` = round(result_item$p_value %||% NA_real_, 4),
        Bandwidth = round(result_item$bandwidth %||% NA_real_, 3),
        `N Left` = result_item$n_left %||% NA_integer_,
        `N Right` = result_item$n_right %||% NA_integer_,
        Error = shorten_strings(result_item$error %||% NA_character_, 40),
        check.names = FALSE
      )
    } else {
      data.frame(`Test Name` = test_name_key, `Matrix Type` = matrix_type_parsed, `Metric Space` = metric_space_parsed,
                 `Test Statistic`=NA, `P-value`=NA, Bandwidth=NA, `N Left`=NA, `N Right`=NA, Error=shorten_strings(as.character(result_item),40), check.names=FALSE)
    }
  }))
}

create_lagged_summary <- function(lagged_io_results_list) {
  if (length(lagged_io_results_list) == 0) {
    return(data.frame())
  }
  summary_rows_list <- list()
  for (result_key_name in names(lagged_io_results_list)) {
    lag_data_item <- lagged_io_results_list[[result_key_name]]
    # RDD results for this lag
    if (!is.null(lag_data_item$rdd_results) && length(lag_data_item$rdd_results) > 0) {
      for (outcome_name_key in names(lag_data_item$rdd_results)) {
        res_item <- lag_data_item$rdd_results[[outcome_name_key]]
        row_id <- paste(result_key_name, "RDD", outcome_name_key, sep = "_")
        if (is.list(res_item)) {
          summary_rows_list[[row_id]] <- data.frame(
            AnalysisID = result_key_name,
            RV_Year = lag_data_item$rv_year %||% NA_integer_,
            IO_Lag = lag_data_item$io_lag_applied %||% NA_integer_,
            IO_Year = lag_data_item$io_year_used %||% NA_integer_,
            ResultType = "RDD Scalar",
            Outcome = stringr::str_to_title(stringr::str_replace_all(res_item$outcome %||% outcome_name_key, "_", " ")),
            Estimate = round(res_item$coefficient %||% NA_real_, 4),
            Std.Error          = round(res_item$se          %||% NA_real_, 4),
            `P-value` = round(res_item$pvalue %||% NA_real_, 4),
            `N Obs (Analysis)` = lag_data_item$n_obs_after_io_processing %||% NA_integer_,
            `N Left (RDD)` = res_item$n_left %||% NA_integer_,
            `N Right (RDD)` = res_item$n_right %||% NA_integer_,
            Error = shorten_strings(res_item$error %||% NA_character_, 30),
            check.names = FALSE
          )
        }
      }
    }
    # Fréchet results for this lag
    if (!is.null(lag_data_item$frechet_results) && length(lag_data_item$frechet_results) > 0) {
      for (test_name_key in names(lag_data_item$frechet_results)) {
        res_item <- lag_data_item$frechet_results[[test_name_key]]
        row_id <- paste(result_key_name, "Frechet", test_name_key, sep = "_")
        if (is.list(res_item)) {
          parts_lag <- stringr::str_split(test_name_key, "_")[[1]]
          matrix_type_lag <- "Unknown Matrix"
          metric_space_lag <- "Unknown Space"
          if (length(parts_lag) >= 2) {
            matrix_type_lag <- paste0(parts_lag[1], ifelse(length(parts_lag)>2 && parts_lag[2]!="covariance" && parts_lag[2]!="network", paste0("_",parts_lag[2]), ""), " matrix")
            metric_space_lag <- parts_lag[length(parts_lag)] 
            if (metric_space_lag=="laplacian") metric_space_lag <- "network (Laplacian)"
          }
          
          summary_rows_list[[row_id]] <- data.frame(
            AnalysisID = result_key_name,
            RV_Year = lag_data_item$rv_year %||% NA_integer_,
            IO_Lag = lag_data_item$io_lag_applied %||% NA_integer_,
            IO_Year = lag_data_item$io_year_used %||% NA_integer_,
            ResultType = paste0("Frechet (", matrix_type_lag, ", ", metric_space_lag, ")"),
            Outcome = test_name_key, # Full Frechet test name
            Estimate = round(res_item$test_statistic %||% NA_real_, 4), # Test Statistic
            Bandwidth = round(res_item$h_variance_used %||% NA_real_, 3) %||% NA_real_,
            `P-value` = round(res_item$p_value %||% NA_real_, 4),
            `N Obs (Analysis)` = lag_data_item$n_obs_after_io_processing %||% NA_integer_,
            `N Left (Frechet)` = res_item$n_left %||% NA_integer_,
            `N Right (Frechet)` = res_item$n_right %||% NA_integer_,
            Error = shorten_strings(res_item$error %||% NA_character_, 30),
            check.names = FALSE
          )
        }
      }
    }
  }
  if (length(summary_rows_list) > 0) {
    return(dplyr::bind_rows(summary_rows_list))
  } else {
    return(data.frame()) # Return empty data.frame if no summary rows generated
  }
}

# Helper to shorten strings for tables
shorten_strings <- function(text, len = 30) {
  if (is.na(text) || text == "" || !nzchar(text)) { # check for empty string too
    return(NA_character_)
  }
  text_char <- as.character(text) # Ensure it's character
  if (nchar(text_char) > len) {
    return(paste0(substr(text_char, 1, len - 3), "..."))
  }
  return(text_char)
}

#### latex table from top industry flows
# Required packages
library(readr)
library(dplyr)
library(knitr)
library(kableExtra)

agg_labels <- c(
  AGR     = "Agri. & Fish.",
  MIN     = "Mining",
  MAN     = "Manufacturing",
  UTL     = "Utilities",
  CON     = "Construction",
  TRD     = "Trade & Hosp.",
  FIN_PRO = "Finance & Business",
  PUB_OTH = "Public & Other"
)

# 1) Read the CSV of top changes
df <- read_csv(file.path(figs, "tables", "leontief_top_changes_rv_2015_io_lag_2.csv"))

# 2) Map the sector codes to full names
df <- df %>%
  mutate(
    Sector_A = agg_labels[Sector_1],
    Sector_B = agg_labels[Sector_2]
  )

# 3) Split into increases and decreases
inc <- df %>% filter(Change_Type == "Increase")
dec <- df %>% filter(Change_Type == "Decrease")

# 4) Build the two subtables, now using Sector_A and Sector_B
inc_tbl <- inc %>%
  select(Rank, Sector_A, Sector_B, Change) %>%
  kable(
    format      = "latex",
    booktabs    = TRUE,
    table.envir = "center",       # emit only the tabular environment
    caption = "Top 5 Increases",
    digits      = 3,
    col.names   = c("Rank", "Sector A", "Sector B", "Change"),
    align       = c("r", "l", "l", "r")
  )

dec_tbl <- dec %>%
  select(Rank, Sector_A, Sector_B, Change) %>%
  kable(
    format      = "latex",
    booktabs    = TRUE,
    table.envir = "center",
    caption = "Top 5 Decreases",
    digits      = 3,
    col.names   = c("Rank", "Sector A", "Sector B", "Change"),
    align       = c("r", "l", "l", "r")
  )

# 5) Assemble full table with two subtables
cat("
\\centering
\\caption{Preferential Tariff Loss: Industry Pairs with Largest Changes}
\\begin{subtable}{0.9\\textwidth}
  \\centering
", inc_tbl, "
\\end{subtable}\\hfill
\\begin{subtable}{0.9\\textwidth}
  \\centering
", dec_tbl, "
\\end{subtable}
", file = file.path(tabs, "leontief_top_changes_rv_2015_io_lag_2.tex"))

# --- RUN THE ANALYSIS ---
# Ensure dataIn and figs are actual paths if not defined elsewhere
# For example:
# dataIn <- "actual/path/to/eora_and_other_data" 
# figs <- "actual/path/to/output_directory"
# config$data_dir <- dataIn
# config$output_dir <- figs
# config$cache_dir <- file.path(figs, "eora_cache")


if ((is.name(config$data_dir) && deparse(config$data_dir) == "dataIn") ||
    (is.name(config$output_dir) && deparse(config$output_dir) == "figs")) {
  warning("`config$data_dir` or `config$output_dir` might be unassigned placeholders (dataIn, figs). Please ensure they are set to actual paths before running.")
  # Example:
  # dataIn <- "/path/to/your/data" # Replace with your actual data path
  # figs <- "/path/to/your/output"  # Replace with your actual output path
  # config$data_dir <- dataIn
  # config$output_dir <- figs
  # config$cache_dir <- file.path(config$output_dir, "eora_cache")
  # cat("Stopping due to placeholder paths. Please define dataIn and figs.\n")
} else {
  if (!dir.exists(config$data_dir)) {
    warning(paste("`config$data_dir` does not exist:", config$data_dir, 
                  "Please ensure the EORA data (e.g., eora_io_data_YYYY folders) is located here."))
  } else {
    cat("Using data_dir:", config$data_dir, "\n")
  }
  if (!dir.exists(config$output_dir)) {
    cat("Output directory will be created at:", config$output_dir, "\n")
    # setup_directories will create it
  } else {
    cat("Using output_dir:", config$output_dir, "\n")
  }
  
  # IMPORTANT: Ensure your frechesTest function is available.
  # If it is from a package (e.g., 'FunctionalRDD'), ensure that package is
  # listed in `required_packages` and can be loaded.
  # The code now tries to use frechesTest::frechesTest or a globally sourced version.
  # If you are using a custom script, you must source it first, e.g.:
  # source("path/to/your_frechesTest_function.R") # Ensure this is done before run_main_analysis
  
  # run_main_analysis(config) # Uncomment to run
}

# Example of how to run if paths are set:
# dataIn <- "C:/Users/username/Documents/eora_project/data_input" # Example path
# figs <- "C:/Users/username/Documents/eora_project/analysis_output" # Example path
config$data_dir <- dataIn
config$output_dir <- figs
config$cache_dir <- file.path(config$output_dir, "eora_cache")
run_main_analysis(config)
