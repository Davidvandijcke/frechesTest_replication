# ================================================================================
# FRÉCHET ANOVA SIMULATIONS
# 
# Purpose: Generate power curves comparing Fréchet Jump Test vs Fréchet ANOVA
#          across different metric spaces (Density, Covariance, Network) and
#          data generating processes (Dynamic, Static)
#
# Author: David Van Dijcke
# ================================================================================

# --- 1. LOAD LIBRARIES -----------------------------------------------------------
library(MASS)
library(Matrix)
library(pracma)
library(osqp)
library(dplyr)
library(knitr)
library(ggplot2)
library(mvtnorm)
library(foreach)
library(doParallel)
library(igraph)
library(tidyr)
library(digest)
library(frechet)
library(RColorBrewer)
library(patchwork)

# Utility operator
`%||%` <- function(a, b) if (is.null(a)) b else a


# --- 2. SIMULATION PARAMETERS -----------------------------------------------------

# Parallel processing settings
parallel_sim <- TRUE  # Set to FALSE for serial debugging

# Simulation parameters
set.seed(123)
N_SIMULATIONS_POWER_CURVE <- 1000  # Number of Monte Carlo simulations per parameter value
ALPHA_LEVEL <- 0.05               # Significance level
CUTOFF_C <- 0.5                   # Jump location
SAMPLE_SIZE_POWER_CURVE <- 200    # Sample size for each simulation

# Parameter grid settings
PARAM_POINTS_ONE_SIDE <- 5
PARAM_SEQ_LENGTH_SYMMETRIC <- 2 * PARAM_POINTS_ONE_SIDE + 1

# Density space parameters: jump in mean
MAX_DENSITY_JUMP <- 2.0
DENSITY_JUMP_PARAM_VALUES <- seq(-MAX_DENSITY_JUMP, MAX_DENSITY_JUMP, 
                                 length.out = PARAM_SEQ_LENGTH_SYMMETRIC)

# Covariance space parameters: multiplicative jump in variance
MAX_LOG_COV_JUMP_FACTOR <- log(2.0)
LOG_COV_JUMP_FACTORS <- seq(-MAX_LOG_COV_JUMP_FACTOR, MAX_LOG_COV_JUMP_FACTOR, 
                            length.out = PARAM_SEQ_LENGTH_SYMMETRIC)
COVARIANCE_JUMP_PARAM_VALUES <- exp(LOG_COV_JUMP_FACTORS)

# Network space parameters: jump in edge probability
MAX_NETWORK_BETA_JUMP <- 0.3
NETWORK_JUMP_PARAM_VALUES <- seq(-MAX_NETWORK_BETA_JUMP, MAX_NETWORK_BETA_JUMP, 
                                 length.out = PARAM_SEQ_LENGTH_SYMMETRIC)

# Data generating process configurations
config_density_dgp <- list(
  BASE_SD = 1.0,              # Standard deviation of observations
  MU_SLOPE_DYNAMIC = 0.8,     # Slope for dynamic DGP
  MU_BASE_STATIC = 0.5        # Base mean for static DGP
)

config_covariance_dgp <- list(
  DIM = 3,                    # Matrix dimension
  SAMPLES_PER_MATRIX = 300,   # Samples to estimate each covariance matrix
  BASE_DIAG_DYNAMIC = 1.5,    # Base diagonal value (dynamic)
  DIAG_SLOPE_DYNAMIC = 0.6,   # Diagonal slope (dynamic)
  BASE_OFFDIAG_DYNAMIC = 0.2, # Base off-diagonal value (dynamic)
  OFFDIAG_SLOPE_DYNAMIC = 0.3,# Off-diagonal slope (dynamic)
  BASE_DIAG_STATIC = 1.5,     # Base diagonal value (static)
  BASE_OFFDIAG_STATIC = 0.2   # Base off-diagonal value (static)
)

config_network_dgp <- list(
  N_NODES = 10,                    # Number of nodes in network
  BETA_CONCENTRATION = 5,          # Beta distribution concentration parameter
  BETA_MEAN_BASE_DYNAMIC = 0.4,    # Base edge probability (dynamic)
  BETA_MEAN_SLOPE_DYNAMIC = 0.2,   # Edge probability slope (dynamic)
  BETA_MEAN_BASE_STATIC = 0.4      # Base edge probability (static)
)


# --- 3. DATA GENERATION FUNCTION --------------------------------------------------

#' Generate data for power curve simulations
#' 
#' @param N_gen Sample size
#' @param X_vals_gen Running variable values
#' @param metric_type_gen One of "Density", "Covariance", "Network"
#' @param dgp_type_gen One of "dynamic", "static"
#' @param current_jump_param_gen Jump magnitude parameter
#' @param cutoff_c_gen Cutoff location
#' @param config_den Configuration for density DGP
#' @param config_cov Configuration for covariance DGP
#' @param config_net Configuration for network DGP
#' @return List with Y_obj (metric space objects) and X_scalar (running variable)
generate_data_for_power_curve <- function(N_gen, X_vals_gen, metric_type_gen, dgp_type_gen,
                                          current_jump_param_gen, cutoff_c_gen,
                                          config_den, config_cov, config_net) {
  
  Y_obj_list <- vector("list", N_gen)
  
  if (metric_type_gen == "Density") {
    # Generate density data: distributions with potential jump in mean
    for (i in 1:N_gen) {
      x_i <- X_vals_gen[i]
      
      # Base mean (no jump)
      if (dgp_type_gen == "dynamic") {
        mu_no_jump <- config_den$MU_SLOPE_DYNAMIC * (x_i - cutoff_c_gen)
      } else {
        mu_no_jump <- config_den$MU_BASE_STATIC
      }
      
      # Add jump if past cutoff
      current_mu <- mu_no_jump
      if (x_i >= cutoff_c_gen) {
        current_mu <- mu_no_jump + current_jump_param_gen
      }
      
      # Generate observations from normal distribution
      Y_obj_list[[i]] <- rnorm(100, mean = current_mu, sd = config_den$BASE_SD)
    }
    
  } else if (metric_type_gen == "Covariance") {
    # Generate covariance matrix data with potential jump in scale
    for (i in 1:N_gen) {
      x_i <- X_vals_gen[i]
      Sigma_base <- diag(config_cov$DIM)
      
      if (dgp_type_gen == "dynamic") {
        # Dynamic DGP: covariance structure changes with X
        diag_val <- config_cov$BASE_DIAG_DYNAMIC + 
          config_cov$DIAG_SLOPE_DYNAMIC * (x_i - cutoff_c_gen)
        offdiag_val <- config_cov$BASE_OFFDIAG_DYNAMIC + 
          config_cov$OFFDIAG_SLOPE_DYNAMIC * (x_i - cutoff_c_gen)
        
        diag(Sigma_base) <- pmax(0.1, diag_val)
        if (config_cov$DIM >= 2) {
          Sigma_base[1,2] <- Sigma_base[2,1] <- offdiag_val
        }
        if (config_cov$DIM >= 3) {
          Sigma_base[1,3] <- Sigma_base[3,1] <- offdiag_val * 0.5
          Sigma_base[2,3] <- Sigma_base[3,2] <- offdiag_val * 0.25
        }
      } else {
        # Static DGP: constant covariance structure
        diag(Sigma_base) <- config_cov$BASE_DIAG_STATIC
        if (config_cov$DIM >= 2) {
          Sigma_base[1,2] <- Sigma_base[2,1] <- config_cov$BASE_OFFDIAG_STATIC
        }
        if (config_cov$DIM >= 3) {
          Sigma_base[1,3] <- Sigma_base[3,1] <- config_cov$BASE_OFFDIAG_STATIC * 0.5
          Sigma_base[2,3] <- Sigma_base[3,2] <- config_cov$BASE_OFFDIAG_STATIC * 0.25
        }
      }
      
      # Ensure positive definiteness
      Sigma_no_jump_pd <- try(
        as.matrix(Matrix::nearPD(Sigma_base, ensureSymmetry = TRUE, 
                                 base.matrix = TRUE)$mat), 
        silent = TRUE
      )
      if(inherits(Sigma_no_jump_pd, "try-error")) {
        Sigma_no_jump_pd <- diag(config_cov$DIM)
      }
      
      # Apply jump if past cutoff
      current_Sigma <- Sigma_no_jump_pd
      if (x_i >= cutoff_c_gen) {
        jump_matrix <- diag(sqrt(current_jump_param_gen), nrow = config_cov$DIM)
        current_Sigma_jumped <- jump_matrix %*% Sigma_no_jump_pd %*% jump_matrix
        
        current_Sigma_jumped_pd <- try(
          as.matrix(Matrix::nearPD(current_Sigma_jumped, ensureSymmetry = TRUE, 
                                   base.matrix = TRUE)$mat), 
          silent = TRUE
        )
        if(!inherits(current_Sigma_jumped_pd, "try-error")) {
          current_Sigma <- current_Sigma_jumped_pd
        }
      }
      
      # Generate data and compute sample covariance
      sample_data_cov <- MASS::mvrnorm(n = config_cov$SAMPLES_PER_MATRIX, 
                                       mu = rep(0, config_cov$DIM), 
                                       Sigma = current_Sigma)
      Y_obj_list[[i]] <- stats::cov(sample_data_cov)
    }
    
  } else if (metric_type_gen == "Network") {
    # Generate network Laplacian data with potential jump in edge probability
    m_nodes <- config_net$N_NODES
    
    for (i in 1:N_gen) {
      x_i <- X_vals_gen[i]
      
      # Base edge probability (no jump)
      if (dgp_type_gen == "dynamic") {
        p_no_jump <- config_net$BETA_MEAN_BASE_DYNAMIC + 
          config_net$BETA_MEAN_SLOPE_DYNAMIC * (x_i - cutoff_c_gen)
      } else {
        p_no_jump <- config_net$BETA_MEAN_BASE_STATIC
      }
      p_no_jump <- pmax(0.05, pmin(0.95, p_no_jump))
      
      # Apply jump if past cutoff
      current_target_p <- p_no_jump
      if (x_i >= cutoff_c_gen) {
        current_target_p <- p_no_jump + current_jump_param_gen
        current_target_p <- pmax(0.05, pmin(0.95, current_target_p))
      }
      
      # Beta distribution parameters for edge weights
      beta_shape1 <- current_target_p * config_net$BETA_CONCENTRATION
      beta_shape2 <- (1 - current_target_p) * config_net$BETA_CONCENTRATION
      if(beta_shape1 <= 0 || beta_shape2 <= 0) {
        beta_shape1 <- 1
        beta_shape2 <- 1
      }
      
      # Generate Laplacian matrix
      L_matrix <- matrix(0, nrow = m_nodes, ncol = m_nodes)
      num_off_diag_upper <- m_nodes * (m_nodes - 1) / 2
      
      if (num_off_diag_upper > 0) {
        beta_samples <- rbeta(num_off_diag_upper, 
                              shape1 = beta_shape1, 
                              shape2 = beta_shape2)
        k_sample <- 1
        for (r_node in 1:(m_nodes - 1)) {
          for (c_node in (r_node + 1):m_nodes) {
            wij <- beta_samples[k_sample]
            L_matrix[r_node, c_node] <- -wij
            L_matrix[c_node, r_node] <- -wij
            k_sample <- k_sample + 1
          }
        }
      }
      
      # Set diagonal to ensure row sums are zero (Laplacian property)
      diag(L_matrix) <- -rowSums(L_matrix)
      Y_obj_list[[i]] <- L_matrix
    }
    
  } else {
    stop(paste("Unsupported metric_type_gen:", metric_type_gen))
  }
  
  return(list(Y_obj = Y_obj_list, X_scalar = X_vals_gen))
}


# --- 4. SINGLE ITERATION FUNCTION -------------------------------------------------

#' Run one power curve iteration
#' 
#' Generates data and applies both Fréchet Jump Test and Fréchet ANOVA
#' 
#' @return Named vector with rejection indicators for both tests
run_one_power_iteration <- function(N_iter, metric_type_iter, dgp_type_iter, 
                                    current_jump_val_iter, cutoff_c_iter, 
                                    alpha_level_iter, config_den_iter, 
                                    config_cov_iter, config_net_iter,
                                    frechesTest_cv_K_folds, 
                                    frechesTest_cv_n_bw_candidates,
                                    frechesTest_verbose) {
  
  # Generate data
  X_vals_iter <- runif(N_iter, 0, 1)
  data_sim_iter <- generate_data_for_power_curve(
    N_iter, X_vals_iter, metric_type_iter, dgp_type_iter,
    current_jump_val_iter, cutoff_c_iter,
    config_den_iter, config_cov_iter, config_net_iter
  )
  Y_obj_sim_iter <- data_sim_iter$Y_obj
  X_scalar_sim_iter <- data_sim_iter$X_scalar
  
  # Set Fréchet options based on metric space
  frechet_opts_iter <- list()
  if (metric_type_iter == "Density") {
    qSup_sim_iter <- seq(0, 1, length.out = 50)
    frechet_opts_iter <- list(
      qSup = qSup_sim_iter,
      den_opts_for_create_density = list(kernel = "gauss", nRegGrid = 50)
    )
  } else if (metric_type_iter == "Covariance") {
    frechet_opts_iter <- list(metric = "frobenius")
  } else if (metric_type_iter == "Network") {
    frechet_opts_iter <- list(metric = "frobenius", W_laplacian_bound = 1.0)
  }
  
  # Apply Fréchet Jump Test
  res_my_test_iter <- tryCatch({
    frechesTest(
      Y_obj = Y_obj_sim_iter, 
      X_scalar = X_scalar_sim_iter, 
      c_val = cutoff_c_iter,
      metric_space_type = tolower(metric_type_iter), 
      h_frechet = "CV", 
      frechet_optns = frechet_opts_iter,
      cv_K_folds = frechesTest_cv_K_folds, 
      cv_n_bw_candidates = frechesTest_cv_n_bw_candidates,
      verbose = frechesTest_verbose
    )
  }, error = function(e) list(p_value = NA_real_, h_variance_used = NA_real_))
  
  pval_my_test_iter <- res_my_test_iter$p_value
  h_selected_iter <- res_my_test_iter$h_variance_used
  
  # Fallback bandwidth if needed
  if (is.na(h_selected_iter) || h_selected_iter <= 1e-5) {
    sd_X_iter <- sd(X_scalar_sim_iter, na.rm = TRUE)
    h_selected_iter <- if (!is.na(sd_X_iter) && sd_X_iter > 1e-6) {
      0.1 * sd_X_iter * (N_iter^(-1/5))
    } else {
      0.05
    }
    h_selected_iter <- max(1e-4, h_selected_iter)
  }
  
  # Apply Fréchet ANOVA (DM test)
  pval_dm_test_iter <- NA_real_
  tryCatch({
    # Select observations within bandwidth of cutoff
    idx_left_dm <- which(X_scalar_sim_iter >= (cutoff_c_iter - h_selected_iter) & 
                           X_scalar_sim_iter < cutoff_c_iter)
    idx_right_dm <- which(X_scalar_sim_iter >= cutoff_c_iter & 
                            X_scalar_sim_iter < (cutoff_c_iter + h_selected_iter))
    
    min_group_size_dm <- 3
    if (length(idx_left_dm) >= min_group_size_dm && 
        length(idx_right_dm) >= min_group_size_dm) {
      
      Y_local_left_dm <- Y_obj_sim_iter[idx_left_dm]
      Y_local_right_dm <- Y_obj_sim_iter[idx_right_dm]
      Y_dm_iter <- c(Y_local_left_dm, Y_local_right_dm)
      group_dm_iter <- c(rep(1, length(Y_local_left_dm)), 
                         rep(2, length(Y_local_right_dm)))
      
      dm_optns_iter_base <- list(boot = FALSE)
      res_dm_iter <- NULL
      
      if (metric_type_iter == "Density") {
        # Special handling for density space
        qSup_for_dm <- frechet_opts_iter$qSup %||% seq(0, 1, length.out = 50)
        den_opts_for_dm_create <- frechet_opts_iter$den_opts_for_create_density %||% 
          list(kernel = "gauss", nRegGrid = 50)
        
        qin_for_dm_test <- tryCatch({
          lapply(Y_dm_iter, function(obs_raw) {
            .get_quantiles_for_obs_jump_test(obs_raw, 
                                             qSup_target = qSup_for_dm, 
                                             den_opts = den_opts_for_dm_create)
          })
        }, error = function(e_q) { NULL })
        
        if (!is.null(qin_for_dm_test) && all(sapply(qin_for_dm_test, is.numeric))) {
          final_dm_opts_den <- dm_optns_iter_base
          res_dm_iter <- frechet::DenANOVA(qin = qin_for_dm_test, 
                                           supin = qSup_for_dm, 
                                           group = group_dm_iter, 
                                           optns = final_dm_opts_den)
        }
        
      } else if (metric_type_iter == "Covariance" || metric_type_iter == "Network") {
        res_dm_iter <- frechet::NetANOVA(Ly = Y_dm_iter, 
                                         group = group_dm_iter, 
                                         optns = dm_optns_iter_base)
      }
      
      if (!is.null(res_dm_iter) && !is.null(res_dm_iter$pvalAsy)) {
        pval_dm_test_iter <- res_dm_iter$pvalAsy
      }
    }
  }, error = function(e) { 
    pval_dm_test_iter <- NA_real_ 
  })
  
  # Return rejection indicators
  return(c(
    my_test_reject = as.integer(!is.na(pval_my_test_iter) && 
                                  pval_my_test_iter < alpha_level_iter),
    dm_test_reject = as.integer(!is.na(pval_dm_test_iter) && 
                                  pval_dm_test_iter < alpha_level_iter)
  ))
}


# --- 5. SETUP PARALLEL BACKEND ----------------------------------------------------

if (parallel_sim) {
  num_cores <- parallel::detectCores() - 1
  if (is.na(num_cores) || num_cores < 1) num_cores <- 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Export utility operator to workers
  clusterEvalQ(cl, {
    `%||%` <- function(a, b) if (is.null(a)) b else a
  })
  
  cat(paste("Starting power curve simulations in PARALLEL on", 
            getDoParWorkers(), "workers.\n"))
} else {
  registerDoSEQ()
  cat("Starting power curve simulations SERIALLY.\n")
}


# --- 6. RUN SIMULATIONS -----------------------------------------------------------

# Initialize results data frame
all_power_results_df <- data.frame()

# Simulation configurations
dgp_types_to_run <- c("dynamic", "static")
metric_spaces_to_run_power <- c("Density", "Covariance", "Network")

# Helper functions to export to parallel workers
frechesTest_helpers_to_export <- c(
  ".get_quantiles_for_obs_jump_test", 
  ".calculate_one_sided_locpoly_weights",
  ".estimate_weighted_frechet_mean", 
  ".estimate_one_sided_frechet_quantities", 
  ".set_bw_range_frechet_style",
  "._cv_for_frechet_jump_h", 
  ".project_to_laplacian_space_osqp"
)

functions_to_export_main <- c(
  "generate_data_for_power_curve", 
  "run_one_power_iteration",
  "frechesTest", 
  frechesTest_helpers_to_export
)

vars_to_export <- c(
  "CUTOFF_C", "ALPHA_LEVEL", "SAMPLE_SIZE_POWER_CURVE",
  "config_density_dgp", "config_covariance_dgp", "config_network_dgp",
  "DENSITY_JUMP_PARAM_VALUES", "COVARIANCE_JUMP_PARAM_VALUES", 
  "NETWORK_JUMP_PARAM_VALUES", "LOG_COV_JUMP_FACTORS"
)

# Run simulations for each combination of DGP type and metric space
for (dgp_type_scenario_loop in dgp_types_to_run) {
  for (metric_scenario_loop in metric_spaces_to_run_power) {
    
    cat(paste("\n--- Generating Power Curve for DGP:", dgp_type_scenario_loop, 
              ", Metric:", metric_scenario_loop, "---\n"))
    
    # Set parameter values and labels based on metric space
    if (metric_scenario_loop == "Density") {
      varying_param_values_current_loop <- DENSITY_JUMP_PARAM_VALUES
      param_name_x_axis_loop <- "Mean Jump (δ)"
    } else if (metric_scenario_loop == "Covariance") {
      varying_param_values_current_loop <- COVARIANCE_JUMP_PARAM_VALUES
      param_name_x_axis_loop <- "Variance Scale Factor (β)"
    } else if (metric_scenario_loop == "Network") {
      varying_param_values_current_loop <- NETWORK_JUMP_PARAM_VALUES
      param_name_x_axis_loop <- "Mean Edge Prob. Jump (Δp)"
    }
    
    # Run simulations (parallel or serial)
    if (parallel_sim) {
      power_results_list_current_scenario <- foreach(
        k_param_idx = seq_along(varying_param_values_current_loop), 
        .combine = 'rbind',
        .packages = c("stats", "MASS", "Matrix", "pracma", "osqp", 
                      "igraph", "frechet", "digest"),
        .export = c(functions_to_export_main, vars_to_export)
      ) %dopar% {
        current_jump_par_val_worker <- varying_param_values_current_loop[k_param_idx]
        
        # Run replications
        rejections_for_param <- t(replicate(N_SIMULATIONS_POWER_CURVE, {
          run_one_power_iteration(
            N_iter = SAMPLE_SIZE_POWER_CURVE, 
            metric_type_iter = metric_scenario_loop, 
            dgp_type_iter = dgp_type_scenario_loop,
            current_jump_val_iter = current_jump_par_val_worker, 
            cutoff_c_iter = CUTOFF_C, 
            alpha_level_iter = ALPHA_LEVEL,
            config_den_iter = config_density_dgp, 
            config_cov_iter = config_covariance_dgp, 
            config_net_iter = config_network_dgp,
            frechesTest_cv_K_folds = 5, 
            frechesTest_cv_n_bw_candidates = 10, 
            frechesTest_verbose = FALSE
          )
        }))
        
        # Calculate power
        data.frame(
          param_value = current_jump_par_val_worker,
          power_my_test = mean(rejections_for_param[, "my_test_reject"], na.rm = TRUE),
          power_dm_test = mean(rejections_for_param[, "dm_test_reject"], na.rm = TRUE)
        )
      }
    } else {
      # Serial execution
      temp_results_serial <- vector("list", length(varying_param_values_current_loop))
      
      for (k_param_idx in seq_along(varying_param_values_current_loop)) {
        cat(paste("  Serial iteration for param value:", 
                  signif(varying_param_values_current_loop[k_param_idx], 3), "\n"))
        
        current_jump_par_val_worker <- varying_param_values_current_loop[k_param_idx]
        
        # Run replications
        rejections_for_param <- t(replicate(N_SIMULATIONS_POWER_CURVE, {
          run_one_power_iteration(
            N_iter = SAMPLE_SIZE_POWER_CURVE, 
            metric_type_iter = metric_scenario_loop, 
            dgp_type_iter = dgp_type_scenario_loop,
            current_jump_val_iter = current_jump_par_val_worker, 
            cutoff_c_iter = CUTOFF_C, 
            alpha_level_iter = ALPHA_LEVEL,
            config_den_iter = config_density_dgp, 
            config_cov_iter = config_covariance_dgp, 
            config_net_iter = config_network_dgp,
            frechesTest_cv_K_folds = 5, 
            frechesTest_cv_n_bw_candidates = 10, 
            frechesTest_verbose = FALSE
          )
        }))
        
        # Calculate power
        temp_results_serial[[k_param_idx]] <- data.frame(
          param_value = current_jump_par_val_worker,
          power_my_test = mean(rejections_for_param[, "my_test_reject"], na.rm = TRUE),
          power_dm_test = mean(rejections_for_param[, "dm_test_reject"], na.rm = TRUE)
        )
      }
      
      power_results_list_current_scenario <- do.call(rbind, temp_results_serial)
    }
    
    # Add metadata to results
    power_results_df_current_scenario <- as.data.frame(power_results_list_current_scenario)
    power_results_df_current_scenario$dgp_type <- dgp_type_scenario_loop
    power_results_df_current_scenario$metric_space <- metric_scenario_loop
    power_results_df_current_scenario$param_name_x_axis <- param_name_x_axis_loop
    
    # Append to overall results
    all_power_results_df <- rbind(all_power_results_df, power_results_df_current_scenario)
  }
}

# Clean up parallel backend
if (parallel_sim) stopCluster(cl)
cat("All power curve simulations complete.\n")


# --- 7. SAVE RESULTS FOR FUTURE USE ----------------------------------------------

# Save the power results dataframe
cat("\n--- Saving simulation results ---\n")
if (!dir.exists("figs/raw_results")) dir.create("figs/raw_results", recursive = TRUE)

saveRDS(all_power_results_df, 
        file = "figs/raw_results/frechetANOVA_power_results.rds")

# Also save key parameters used in the simulation
simulation_params <- list(
  N_SIMULATIONS_POWER_CURVE = N_SIMULATIONS_POWER_CURVE,
  ALPHA_LEVEL = ALPHA_LEVEL,
  CUTOFF_C = CUTOFF_C,
  SAMPLE_SIZE_POWER_CURVE = SAMPLE_SIZE_POWER_CURVE,
  DENSITY_JUMP_PARAM_VALUES = DENSITY_JUMP_PARAM_VALUES,
  COVARIANCE_JUMP_PARAM_VALUES = COVARIANCE_JUMP_PARAM_VALUES,
  NETWORK_JUMP_PARAM_VALUES = NETWORK_JUMP_PARAM_VALUES,
  LOG_COV_JUMP_FACTORS = LOG_COV_JUMP_FACTORS,
  dgp_types_to_run = dgp_types_to_run,
  metric_spaces_to_run_power = metric_spaces_to_run_power
)

saveRDS(simulation_params, 
        file = "figs/raw_results/frechetANOVA_simulation_params.rds")

cat("Power results and parameters saved to figs/raw_results/\n")


# --- 8. GENERATE TRANSPOSED POWER CURVES PLOT (2x3 instead of 3x2) ---------------

if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("tidyr", quietly = TRUE) &&
    requireNamespace("patchwork", quietly = TRUE)) {
  
  cat("\n--- Generating Combined Power Curve Plot (2x3 Layout) ---\n")
  
  # Initialize plot list
  plot_list_final <- list()
  subplot_idx_counter <- 0
  
  # Define test aesthetics
  test_colors <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1:2]
  test_shapes <- c(16, 17)
  test_linetypes <- c("solid", "longdash")
  test_labels <- c("Fréchet Jump Test", "Fréchet ANOVA")
  
  names(test_colors) <- test_labels
  names(test_shapes) <- test_labels
  names(test_linetypes) <- test_labels
  
  # Create plots: DGP types in rows, metric spaces in columns
  for (dgp_plot_loop in dgp_types_to_run) {
    for (metric_plot_loop in metric_spaces_to_run_power) {
      subplot_idx_counter <- subplot_idx_counter + 1
      
      # Filter data for current subplot
      current_plot_data_df_plot <- all_power_results_df %>%
        filter(metric_space == metric_plot_loop, dgp_type == dgp_plot_loop)
      
      if(nrow(current_plot_data_df_plot) == 0) {
        plot_list_final[[subplot_idx_counter]] <- patchwork::plot_spacer()
        next
      }
      
      # Extract x-axis label
      current_x_axis_label_plot <- unique(current_plot_data_df_plot$param_name_x_axis)
      
      # Create subplot title
      dgp_label <- ifelse(dgp_plot_loop == "dynamic", "Pw-Smooth", "Pw-Const.")
      subplot_title_text_full <- paste0("(", LETTERS[subplot_idx_counter], ") ",
                                        tools::toTitleCase(metric_plot_loop),
                                        " (", dgp_label, " DGP)")
      # Reshape data for plotting
      power_results_long_plot <- tidyr::pivot_longer(
        current_plot_data_df_plot, 
        cols = c("power_my_test", "power_dm_test"), 
        names_to = "Test_raw", 
        values_to = "Power"
      )
      
      power_results_long_plot$Test <- factor(
        power_results_long_plot$Test_raw, 
        levels = c("power_my_test", "power_dm_test"), 
        labels = test_labels
      )
      
      # Create subplot with LARGER text sizes
      p_current_subplot <- ggplot(power_results_long_plot, 
                                  aes(x = param_value, y = Power, 
                                      color = Test, linetype = Test, shape = Test)) +
        geom_line(linewidth = 0.9) + 
        geom_point(size = 2.5, stroke = 0.7) +
        geom_hline(yintercept = ALPHA_LEVEL, linetype = "dashed", 
                   color = "grey30", linewidth = 0.7) +
        labs(title = subplot_title_text_full, 
             x = current_x_axis_label_plot, 
             y = NULL) +
        scale_y_continuous(limits = c(-0.02, 1.02), 
                           breaks = seq(0, 1, 0.2), 
                           expand = c(0, 0.015)) +
        scale_color_manual(values = test_colors, name = "Test Method") + 
        scale_shape_manual(values = test_shapes, name = "Test Method") +
        scale_linetype_manual(values = test_linetypes, name = "Test Method") +
        theme_classic(base_size = 14) +  # Increased from 11 to 14
        theme(
          legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = rel(1.0), face = "plain"),  # Larger title
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"), 
          axis.title.x = element_text(size = rel(0.9), margin = margin(t = 5)),  # Larger x-axis label
          axis.title.y = element_text(size = rel(0.9), margin = margin(r = 5)),  # Larger y-axis label
          axis.text = element_text(size = rel(0.8))  # Larger axis text
        )
      
      # Special x-axis handling for covariance space (log scale)
      if (metric_plot_loop == "Covariance") {
        log_breaks_plot_subplot <- pretty(LOG_COV_JUMP_FACTORS, n = 4)
        actual_breaks_plot_subplot <- exp(log_breaks_plot_subplot)
        actual_labels_plot_subplot <- sapply(actual_breaks_plot_subplot, function(br) {
          if (abs(log(br)) < 1e-6) "1.0" else sprintf("%.1f", br)
        })
        p_current_subplot <- p_current_subplot + 
          scale_x_continuous(breaks = actual_breaks_plot_subplot, 
                             labels = actual_labels_plot_subplot)
      } else {
        p_current_subplot <- p_current_subplot + 
          scale_x_continuous(n.breaks = 5)
      }
      
      # Add y-axis label only to first column
      if(metric_plot_loop == metric_spaces_to_run_power[1]) {
        p_current_subplot <- p_current_subplot + 
          labs(y = "Empirical Power")
      } else {
        p_current_subplot <- p_current_subplot + 
          theme(axis.text.y = element_blank(), 
                axis.ticks.y = element_blank(), 
                axis.title.y = element_blank())
      }
      
      plot_list_final[[subplot_idx_counter]] <- p_current_subplot
    }
  }
  
  # Combine plots in 2x3 layout (DGP types in rows, metric spaces in columns)
  if (length(plot_list_final) == 6) {
    combined_figure_patch <- (plot_list_final[[1]] | plot_list_final[[2]] | plot_list_final[[3]]) / 
      (plot_list_final[[4]] | plot_list_final[[5]] | plot_list_final[[6]])
    
    # Add legend and caption with CLEAN legend box
    final_figure_with_legend <- combined_figure_patch + 
      plot_layout(guides = "collect")  &
      theme(
        legend.position = "bottom", 
        legend.box.margin = margin(t = 15, b = 5), 
        legend.title = element_text(face = "bold", size = rel(1.1)),  # Larger legend title
        legend.text = element_text(size = rel(1.0)),  # Larger legend text
        plot.caption = element_text(hjust = 0.5, size = rel(1.0),  # Larger caption
                                    margin = margin(t = 15, b = 5)),
        # Remove grey box around legend
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)
      )
    
    print(final_figure_with_legend)
    
    # Save plot
    figs_dir_path <- figs 
    if (!dir.exists(figs_dir_path)) dir.create(figs_dir_path, recursive = TRUE)
    
    combined_plot_filename_final <- paste0("power_curves_stacked_2x3_N", 
                                           SAMPLE_SIZE_POWER_CURVE, ".png")
    
    ggsave(file.path(figs_dir_path, combined_plot_filename_final), 
           plot = final_figure_with_legend, 
           width = 12, height = 8, dpi = 300,
           bg = "white")  # Ensure white background in saved file
    
    cat(paste0("Combined stacked power curve plot (2x3 layout) saved to ", 
               file.path(figs_dir_path, combined_plot_filename_final), "\n"))
  } else {
    cat("Incorrect number of plots generated for stacking.\n")
  }
} else {
  cat("\nPackages ggplot2, tidyr, and patchwork not found.\n")
  print(knitr::kable(all_power_results_df, digits = 3, format = "pipe", 
                     caption = paste("Full Power Comparison Results at N =", 
                                     SAMPLE_SIZE_POWER_CURVE)))
}