# Needed libraries (ensure they are installed)
library(MASS)
library(Matrix)
library(pracma)
library(osqp)
library(dplyr)
library(knitr)
# library(ggplot2) # Not strictly needed for this script if only generating tables
library(mvtnorm)
library(foreach)
library(doParallel)
library(igraph)
library(xtable)
library(data.table)
library(tidyr) # For pivot_wider in table prep

# --- Source Your FrechesTest Code ---
`%||%` <- function(a, b) if (is.null(a)) b else a
# source("FrechetJumpTest_helpers.R") # Ensure these are available
# source("frechesTest.R")

# --- Simulation Parameters ---
parallel_main_sim <- TRUE # Set to FALSE for serial debugging
metric_spaces_to_run <- c("Density", "Covariance", "Network")
N_SIMULATIONS <- 1000 # Number of Monte Carlo repetitions for the paper
ALPHA_LEVEL <- 0.05
CUTOFF_C <- 0.5
SAMPLE_SIZES <- c(200, 500, 1000)

# --- DGP Parameters for H0 (Smooth component) and H1 (Jump Magnitude) ---
# These are now more clearly separated for use in the generate_* functions

# Density Space
DGP_PARAMS_DENSITY <- list(
  BASE_SD = 1.0,
  MU_SLOPE = 0.8,          # Slope for the smooth part E[Y|X=x] = MU_SLOPE * (x - CUTOFF_C)
  MU_JUMP_H1 = 1.5         # Additive jump for H1
)

# Covariance Space
DGP_PARAMS_COVARIANCE <- list(
  DIM = 3,
  SAMPLES_PER_MATRIX = 300,
  BASE_DIAG_SLOPE = 0.6,   # Slope for diagonal elements
  BASE_DIAG_INTERCEPT = 1.5, # Intercept for diagonal at x = CUTOFF_C
  OFFDIAG_SLOPE = 0.3,     # Slope for primary off-diagonal
  OFFDIAG_INTERCEPT = 0.2, # Intercept for primary off-diagonal at x = CUTOFF_C
  JUMP_FACTOR_H1 = 1.5     # Multiplicative scaling factor for SDs under H1
)

# Network Space (ZM-style Beta-Laplacian)
DGP_PARAMS_NETWORK <- list(
  N_NODES = 10,
  BETA_MEAN_SLOPE = 0.2,      # Slope for mean edge probability p(x)
  BETA_MEAN_INTERCEPT = 0.4,  # Intercept for mean p(x) at x = CUTOFF_C
  BETA_JUMP_H1 = 0.25,        # Additive jump to mean p(x) for H1
  BETA_CONCENTRATION = 5
)


# --- Helper Functions for Data Generation (Harmonized Structure) ---
generate_density_data_main <- function(N, X_vals, has_jump_flag,
                                       cutoff_c_val = CUTOFF_C,
                                       dgp_p = DGP_PARAMS_DENSITY) {
  Y_obj_list <- vector("list", N)
  for (i in 1:N) {
    x_i <- X_vals[i]
    # Smooth component: E[Y|X=x] = slope * (x - c) (intercept is implicitly 0 at c)
    # If you want a non-zero intercept for the smooth part at x=c, add it here.
    # For simplicity, let's assume the MU_SLOPE defines change relative to value at c=0.
    # Or, to match your power curve: smooth part is slope * (x-c) + base_static_mu_at_c_for_dynamic
    # Let's make it simpler: the slope defines the change, and any "base" value at c is 0 for the smooth part.
    # The description says: mu(x) = 0.8(x-0.5). This implies intercept at x=0.5 is 0 for the smooth part.
    mu_smooth <- dgp_p$MU_SLOPE * (x_i - cutoff_c_val) # Dynamic component
    
    current_mu <- mu_smooth
    if (has_jump_flag && x_i >= cutoff_c_val) {
      current_mu <- mu_smooth + dgp_p$MU_JUMP_H1 # Add fixed H1 jump
    }
    Y_obj_list[[i]] <- rnorm(100, mean = current_mu, sd = dgp_p$BASE_SD)
  }
  return(list(Y_obj = Y_obj_list, X_scalar = X_vals))
}

generate_covariance_data_main <- function(N, X_vals, has_jump_flag,
                                          cutoff_c_val = CUTOFF_C,
                                          dgp_p = DGP_PARAMS_COVARIANCE) {
  Y_obj_list <- vector("list", N)
  for (i in 1:N) {
    x_i <- X_vals[i]
    # Smooth component for Sigma(x)
    diag_val_smooth <- dgp_p$BASE_DIAG_INTERCEPT + dgp_p$BASE_DIAG_SLOPE * (x_i - cutoff_c_val)
    offdiag_val_smooth <- dgp_p$OFFDIAG_INTERCEPT + dgp_p$OFFDIAG_SLOPE * (x_i - cutoff_c_val)
    
    Sigma_base <- diag(dgp_p$DIM)
    diag(Sigma_base) <- pmax(0.1, diag_val_smooth)
    if(dgp_p$DIM >= 2) Sigma_base[1,2] <- Sigma_base[2,1] <- offdiag_val_smooth
    if(dgp_p$DIM >= 3) Sigma_base[1,3] <- Sigma_base[3,1] <- offdiag_val_smooth * 0.5
    if(dgp_p$DIM >= 3) Sigma_base[2,3] <- Sigma_base[3,2] <- offdiag_val_smooth * 0.25
    
    Sigma_no_jump_pd <- try(as.matrix(Matrix::nearPD(Sigma_base, ensureSymmetry = TRUE, base.matrix=TRUE)$mat), silent=TRUE)
    if(inherits(Sigma_no_jump_pd, "try-error")) Sigma_no_jump_pd <- diag(dgp_p$DIM)
    
    current_Sigma <- Sigma_no_jump_pd
    if (has_jump_flag && x_i >= cutoff_c_val) {
      jump_diag_scale_matrix <- diag(sqrt(dgp_p$JUMP_FACTOR_H1), nrow = dgp_p$DIM)
      current_Sigma_jumped <- jump_diag_scale_matrix %*% Sigma_no_jump_pd %*% jump_diag_scale_matrix
      current_Sigma_jumped_pd <- try(as.matrix(Matrix::nearPD(current_Sigma_jumped, ensureSymmetry = TRUE, base.matrix=TRUE)$mat), silent=TRUE)
      if(!inherits(current_Sigma_jumped_pd, "try-error")) current_Sigma <- current_Sigma_jumped_pd
    }
    sample_data_cov <- MASS::mvrnorm(n = dgp_p$SAMPLES_PER_MATRIX, mu = rep(0, dgp_p$DIM), Sigma = current_Sigma)
    Y_obj_list[[i]] <- stats::cov(sample_data_cov)
  }
  return(list(Y_obj = Y_obj_list, X_scalar = X_vals))
}

generate_network_data_main <- function(N, X_vals, has_jump_flag,
                                       cutoff_c_val = CUTOFF_C,
                                       dgp_p = DGP_PARAMS_NETWORK) {
  Y_obj_list <- vector("list", N)
  for (i in 1:N) {
    x_i <- X_vals[i]
    # Smooth component for mean edge probability p(x)
    p_smooth <- dgp_p$BETA_MEAN_INTERCEPT + dgp_p$BETA_MEAN_SLOPE * (x_i - cutoff_c_val)
    p_smooth <- pmax(0.05, pmin(0.95, p_smooth))
    
    current_target_p <- p_smooth
    if (has_jump_flag && x_i >= cutoff_c_val) {
      current_target_p <- p_smooth + dgp_p$BETA_JUMP_H1 # Add fixed H1 jump
      current_target_p <- pmax(0.05, pmin(0.95, current_target_p))
    }
    beta_shape1 <- current_target_p * dgp_p$BETA_CONCENTRATION
    beta_shape2 <- (1 - current_target_p) * dgp_p$BETA_CONCENTRATION
    if(beta_shape1 <=0 || beta_shape2 <=0) {beta_shape1 <- 1; beta_shape2 <- 1;}
    L_matrix <- matrix(0, nrow = dgp_p$N_NODES, ncol = dgp_p$N_NODES)
    num_off_diag_upper <- dgp_p$N_NODES * (dgp_p$N_NODES - 1) / 2
    if (num_off_diag_upper > 0) {
      beta_samples <- rbeta(num_off_diag_upper, shape1 = beta_shape1, shape2 = beta_shape2)
      k_sample <- 1
      for (r_node in 1:(dgp_p$N_NODES - 1)) for (c_node in (r_node + 1):dgp_p$N_NODES) {
        wij <- beta_samples[k_sample]; L_matrix[r_node, c_node] <- -wij; L_matrix[c_node, r_node] <- -wij; k_sample <- k_sample + 1
      }
    }
    diag(L_matrix) <- -rowSums(L_matrix)
    Y_obj_list[[i]] <- L_matrix
  }
  return(list(Y_obj = Y_obj_list, X_scalar = X_vals))
}


# --- Main Simulation Loop Function (Modified to use new DGP functions) ---
run_single_simulation_set <- function(N_val_iter, metric_val_iter, has_jump_iter_flag) {
  # X_vals generation
  X_vals <- runif(N_val_iter, 0, 1) # For main sim, X can span [0,1] for H0/H1
  # Or, if specific RDD setup needed for H1:
  # if (has_jump_iter_flag) {
  #   n_half <- N_val_iter / 2
  #   X_vals <- c(runif(ceiling(n_half), 0, CUTOFF_C - 1e-6), # ensure strictly less for clarity
  #               runif(floor(n_half), CUTOFF_C, 1))
  # } else {
  #   X_vals <- runif(N_val_iter, 0, 1)
  # }
  
  
  data_list <- NULL
  frechet_opts_config <- list()
  
  if (metric_val_iter == "Density") {
    data_list <- generate_density_data_main(N_val_iter, X_vals, has_jump_iter_flag)
    qSup_sim <- seq(0, 1, length.out = 50) # Consistent with power curve
    frechet_opts_config <- list(qSup = qSup_sim, den_opts_for_create_density = list(kernel = "gauss", nRegGrid = 50))
  } else if (metric_val_iter == "Covariance") {
    data_list <- generate_covariance_data_main(N_val_iter, X_vals, has_jump_iter_flag)
    frechet_opts_config <- list(metric = "frobenius")
  } else if (metric_val_iter == "Network") {
    data_list <- generate_network_data_main(N_val_iter, X_vals, has_jump_iter_flag)
    frechet_opts_config <- list(metric = "frobenius", W_laplacian_bound = 1.0)
  } else {
    stop(paste("Unknown metric space for simulation:", metric_val_iter))
  }
  
  if (is.null(data_list) || is.null(data_list$Y_obj) || length(data_list$Y_obj) == 0) {
    warning(paste("Data generation failed for N=", N_val_iter, "Metric=", metric_val_iter))
    return(NA_real_)
  }
  
  result <- tryCatch({
    frechesTest(
      Y_obj = data_list$Y_obj, X_scalar = data_list$X_scalar, c_val = CUTOFF_C,
      metric_space_type = tolower(metric_val_iter),
      h_frechet = "CV", kernel_frechet_char = "epan",
      frechet_optns = frechet_opts_config,
      cv_K_folds = 5, cv_n_bw_candidates = 10,
      min_bw_cv_factor = 0.001, undersmooth_factor=0.8,
      verbose = FALSE
    )
  }, error = function(e) {
    return(list(p_value = NA_real_))
  })
  return(result$p_value)
}


# --- Setup Parallel Backend (Conditional) ---
if (parallel_main_sim) {
  num_cores <- detectCores() - 1
  if (is.na(num_cores) || num_cores < 1) num_cores <- 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    `%||%` <- function(a, b) if (is.null(a)) b else a
    # Source files if they are not part of a package loaded by .packages
    # source("FrechetJumpTest_helpers.R")
    # source("frechesTest.R")
  })
  cat(paste("Running main simulations in PARALLEL on", getDoParWorkers(), "workers.\n"))
} else {
  registerDoSEQ()
  cat(paste("Running main simulations SERIALLY.\n"))
}


# --- Storing Results ---
results_list <- list()

# --- Iterate Over Scenarios ---
frechesTest_helpers_to_export_main_sim <- c(
  ".get_quantiles_for_obs_jump_test", ".calculate_one_sided_locpoly_weights",
  ".estimate_weighted_frechet_mean", ".estimate_one_sided_frechet_quantities",
  ".set_bw_range_frechet_style", "._cv_for_frechet_jump_h",
  ".project_to_laplacian_space_osqp"
)
functions_to_export_main_sim <- c("run_single_simulation_set", "frechesTest",
                                  "generate_density_data_main", "generate_covariance_data_main", "generate_network_data_main", # Use new _main versions
                                  frechesTest_helpers_to_export_main_sim)
# These are now lists of parameters
dgp_configs_to_export <- c("DGP_PARAMS_DENSITY", "DGP_PARAMS_COVARIANCE", "DGP_PARAMS_NETWORK")
vars_to_export_main_sim <- c("CUTOFF_C", "ALPHA_LEVEL", dgp_configs_to_export)


for (N_val_main in SAMPLE_SIZES) {
  for (metric_val_main in metric_spaces_to_run) {
    cat(paste0("\nScenario: N=", N_val_main, ", Metric=", metric_val_main, "\n"))
    
    # H0: No Jump
    cat("  Running H0 (No Jump)...\n")
    if (parallel_main_sim) {
      p_values_h0 <- foreach(
        i = 1:N_SIMULATIONS, .combine = 'c',
        .packages = c("MASS", "Matrix", "stats", "pracma", "osqp", "igraph", "digest", "frechet"),
        .export = c(functions_to_export_main_sim, vars_to_export_main_sim), # Export config lists
        .errorhandling = 'pass'
      ) %dopar% {
        run_single_simulation_set(N_val_main, metric_val_main, FALSE)
      }
    } else { # Serial execution for H0
      p_values_h0_list <- vector("list", N_SIMULATIONS)
      for(i_sim in 1:N_SIMULATIONS) { # Renamed loop variable
        if(i_sim %% (N_SIMULATIONS/10) == 0) cat(".")
        p_values_h0_list[[i_sim]] <- run_single_simulation_set(N_val_main, metric_val_main, FALSE)
      }
      p_values_h0 <- unlist(p_values_h0_list)
      cat("\n")
    }
    num_errors_h0 <- sum(sapply(p_values_h0, function(x) inherits(x, "error") || inherits(x, "simpleError")))
    p_values_h0_clean <- p_values_h0[!sapply(p_values_h0, function(x) inherits(x, "error") || inherits(x, "simpleError"))]
    p_values_h0_clean <- as.numeric(unlist(p_values_h0_clean))
    rejection_rate_h0 <- mean(p_values_h0_clean < ALPHA_LEVEL, na.rm = TRUE)
    num_na_h0 <- sum(is.na(p_values_h0_clean)) + num_errors_h0 # num_errors already counts errors
    results_list[[length(results_list) + 1]] <- data.frame(
      SampleSize = N_val_main, MetricSpace = metric_val_main, Hypothesis = "H0",
      RejectionRate = rejection_rate_h0, NumNA = num_na_h0 # Removed NumErrors as it's part of NumNA
    )
    
    # H1: With Jump
    cat("  Running H1 (With Jump)...\n")
    if (parallel_main_sim) {
      p_values_h1 <- foreach(
        i = 1:N_SIMULATIONS, .combine = 'c',
        .packages = c("MASS", "Matrix", "stats", "pracma", "osqp", "igraph", "digest", "frechet"),
        .export = c(functions_to_export_main_sim, vars_to_export_main_sim),
        .errorhandling = 'pass'
      ) %dopar% {
        run_single_simulation_set(N_val_main, metric_val_main, TRUE)
      }
    } else { # Serial execution for H1
      p_values_h1_list <- vector("list", N_SIMULATIONS)
      for(i_sim in 1:N_SIMULATIONS) {
        if(i_sim %% (N_SIMULATIONS/10) == 0) cat(".")
        p_values_h1_list[[i_sim]] <- run_single_simulation_set(N_val_main, metric_val_main, TRUE)
      }
      p_values_h1 <- unlist(p_values_h1_list)
      cat("\n")
    }
    num_errors_h1 <- sum(sapply(p_values_h1, function(x) inherits(x, "error") || inherits(x, "simpleError")))
    p_values_h1_clean <- p_values_h1[!sapply(p_values_h1, function(x) inherits(x, "error") || inherits(x, "simpleError"))]
    p_values_h1_clean <- as.numeric(unlist(p_values_h1_clean))
    rejection_rate_h1 <- mean(p_values_h1_clean < ALPHA_LEVEL, na.rm = TRUE)
    num_na_h1 <- sum(is.na(p_values_h1_clean)) + num_errors_h1
    results_list[[length(results_list) + 1]] <- data.frame(
      SampleSize = N_val_main, MetricSpace = metric_val_main, Hypothesis = "H1",
      RejectionRate = rejection_rate_h1, NumNA = num_na_h1
    )
  }
}

# --- Stop Parallel Backend ---
if (parallel_main_sim) stopCluster(cl)

# Combine results
final_results_df <- do.call(rbind, results_list)

fwrite(final_results_df, file=file.path(dataOut, "sims_table_df.csv"))

# --- Display and Generate LaTeX Table ---
cat("\n\n--- Simulation Results Summary ---\n")
cat("Nominal Alpha Level:", ALPHA_LEVEL, "\n")
cat("Number of Simulations per Scenario:", N_SIMULATIONS, "\n\n")

table_for_latex_prep <- final_results_df %>%
  mutate(
    RejectionRate_Display = sprintf("%.3f", round(RejectionRate, 3)),
    NumNA_Display = as.integer(NumNA) # Ensure NumNA is integer for display
  ) %>%
  select(MetricSpace, SampleSize, Hypothesis, RejectionRate_Display, NumNA_Display)

# Pivot for H0 and H1 side-by-side for RejectionRate and NumNA
table_for_latex <- tidyr::pivot_wider(
  table_for_latex_prep,
  names_from = Hypothesis,
  values_from = c(RejectionRate_Display, NumNA_Display),
  names_glue = "{Hypothesis}_{.value}"
) %>%
  select(
    MetricSpace, SampleSize,
    `H0_RejectionRate_Display`, `H1_RejectionRate_Display`
  ) %>%
  arrange(MetricSpace, SampleSize)

colnames(table_for_latex) <- c(
  "Metric Space", "N",
  "Size", "Power",
  "NAs (H$_0$)", "NAs (H$_1$)" # Adjusted column names
)

cat("\n--- LaTeX Table Code ---\n")

addtorow_econ <- NULL # No internal midrules for this style

# Create the xtable object. Caption and label are defined here.
xtable_obj_econ <- xtable::xtable(table_for_latex, # table_for_latex should now have 4 data columns
                                  caption = paste0("Simulation Study: Empirical Size and Power ($\\alpha = ", ALPHA_LEVEL,", ", N_SIMULATIONS, "$ Repetitions per Scenario)"),
                                  label = "tab:simulation_results_econ", # New label
                                  align = "llccr",  # l (rownames), l (Metric), c (N), c (Size), r (Power) - adjust as needed
                                  # For example, to make Size and Power align by decimal with siunitx later: "llS[table-format=1.3]S[table-format=1.3]"
                                  # But for basic xtable, c or r is common.
                                  digits = c(0, 0, 0, 3, 3) # Digits: rownames, Metric, N, Size, Power
)

# Generate the LaTeX string for the tabular environment
latex_tabular_string_econ <- print(xtable_obj_econ,
                                   type = "latex",
                                   floating = FALSE,
                                   include.rownames = FALSE,
                                   booktabs = TRUE,               # ESSENTIAL for \toprule, \midrule, \bottomrule
                                   #hline.after = NULL,            # REMOVE all default hlines from xtable; booktabs handles it
                                   add.to.row = addtorow_econ,    # Set to NULL for strict econ style
                                   sanitize.text.function = function(x){x},
                                   print.results = FALSE
)
latex_tabular_string_econ

# Define the output file path
output_tex_tabular_file_econ <- "simulation_results_tabular.tex"

# Write the string to the file
writeLines(latex_tabular_string_econ, con = file.path(tabs, output_tex_tabular_file_econ))
