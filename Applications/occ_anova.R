################################################################################
##  SIPP-based Fréchet RDD Analysis
##  Multi-state, multi-outcome version
##  Including Oregon 2008-2021 moving pay threshold
################################################################################

# Load panel data
panel <- fread(file=file.path(dataOut, "sipp_job_panel.csv"))

################################################################################
##  0. Analysis Configuration
################################################################################

# Analysis parameters - Fixed to WA state only as in paper
STATES  <- c("WA")         # Washington state only
OUTCOME <- "labor"         # Options: "occ4", "ind4", "labor", "class"
IND_LEVEL <- 1L            # 1=legacy 4 clusters, 2=NAICS 2-digit, 3=3-digit, 4=4-digit
OCC_LEVEL <- 1L            # 1=legacy 4 SOC buckets, 2=SOC 2-digit, 3=3-digit, 4=4-digit

################################################################################
##  1. State Configuration Lookup Table
################################################################################

# State configuration - Washington state only
state_cfg <- list(
  WA = list(
    fips = "53",
    yrs  = 2020:2023,
    wage = "annual_salary",
    cuts = data.table(
      end_year = c(2020, 2021, 2022, 2023),
      cutoff   = c(100000, 101390, 107301.04, 116593.18)
    ),
    scale = 1
  )
)

################################################################################
##  2. Build Pooled Configuration
################################################################################

cfg_list <- state_cfg[STATES]
if (any(vapply(cfg_list, is.null, logical(1))))
  stop("One of the STATE codes in STATES is not in state_cfg.")

cfg <- list(
  fips  = vapply(cfg_list, `[[`, "", "fips"),
  yrs   = unique(unlist(lapply(cfg_list, `[[`, "yrs"))),
  wage  = unique(vapply(cfg_list, `[[`, "",  "wage")),
  scale = unique(vapply(cfg_list, `[[`, 1 , "scale"))
)

if (length(cfg$wage ) > 1)
  stop("Mixed running-variable definitions (hourly vs annual).")
if (length(cfg$scale) > 1)
  stop("Mixed scaling factors – please use comparable states.")

cfg$cuts <- rbindlist(
  lapply(cfg_list, \(l) copy(l$cuts)[, fips := l$fips ])
)

################################################################################
##  3. Filter Panel Data and Merge Cutoffs
################################################################################

panel[, state_now := sprintf("%02d", as.integer(state_now))]

win_start <- min(cfg$yrs)
win_end   <- max(cfg$yrs)

dat <- panel[
  state_now %in% cfg$fips &
    start_year >= win_start &     
    !is.na(get(cfg$wage))
]

dat <- dat[class_worker_now %in% c(5,6)]

# Helper function: identify salaried-exempt workers (SOC 2-digit codes 11-29)
is_salaried_exempt <- function(soc_code){
  soc2 <- suppressWarnings(as.integer(substr(sprintf("%06s", soc_code), 1, 2)))
  !is.na(soc2) & soc2 >= 11 & soc2 <= 29
}

dat[, salaried_exempt := is_salaried_exempt(occupation_now)]

# Merge cutoff values for Washington state
dat <- merge(
  dat,
  cfg$cuts,
  by.x = c("state_now", "end_year"),
  by.y = c("fips", "end_year"),
  all.x = TRUE,
  sort  = FALSE
)

# Create centered running variable
dat[, X := (get(cfg$wage) - cutoff) / cfg$scale]

# Keep clean copy for diagnostics
temp <- copy(dat)

################################################################################
## 4 Generate Summary Statistics Table for LaTeX
################################################################################

# Create summary statistics for key variables
summary_vars <- c("annual_salary", "hourly_wage", "age_start", "educ_level", 
                  "wfh_0", "wfh_1", "wfh_2", "wfh_3", "wfh_4", "wfh_5",
                  "salaried_exempt")

# Function to calculate summary stats
calc_summary_stats <- function(dt, vars) {
  stats_list <- lapply(vars, function(v) {
    if (v %in% names(dt) && sum(!is.na(dt[[v]])) > 0) {
      x <- dt[[v]]
      if (is.numeric(x)) {
        data.table(
          Variable = v,
          N = sum(!is.na(x)),
          Mean = mean(x, na.rm = TRUE),
          SD = sd(x, na.rm = TRUE),
          Min = min(x, na.rm = TRUE),
          P25 = quantile(x, 0.25, na.rm = TRUE),
          Median = median(x, na.rm = TRUE),
          P75 = quantile(x, 0.75, na.rm = TRUE),
          Max = max(x, na.rm = TRUE)
        )
      }
    }
  })
  rbindlist(stats_list[!sapply(stats_list, is.null)])
}

# Calculate summary statistics
summary_stats <- calc_summary_stats(dat, summary_vars)

# Create LaTeX table
create_latex_table <- function(stats_dt) {
  # Variable labels for prettier output
  var_labels <- c(
    annual_salary = "Annual Salary (\\$)",
    hourly_wage = "Hourly Wage (\\$/hr)",
    age_start = "Age",
    educ_level = "Education Level",
    wfh_0 = "WFH 0 days/week",
    wfh_1 = "WFH 1 day/week",
    wfh_2 = "WFH 2 days/week",
    wfh_3 = "WFH 3 days/week",
    wfh_4 = "WFH 4 days/week",
    wfh_5 = "WFH 5 days/week",
    salaried_exempt = "Salaried Exempt"
  )
  
  # Start building LaTeX table
  latex_lines <- c(
    "\\begin{tabular}{lrrrrrrrr}",
    "\\toprule",
    "Variable & N & Mean & SD & Min & P25 & Median & P75 & Max \\\\",
    "\\midrule"
  )
  
  # Add data rows
  for (i in 1:nrow(stats_dt)) {
    var_name <- stats_dt$Variable[i]
    label <- ifelse(var_name %in% names(var_labels), 
                    var_labels[[var_name]], 
                    var_name)
    
    # Format values based on variable type
    if (var_name %in% c("annual_salary")) {
      row_values <- c(
        format_number(stats_dt$N[i], 0),
        format_number(stats_dt$Mean[i], 0),
        format_number(stats_dt$SD[i], 0),
        format_number(stats_dt$Min[i], 0),
        format_number(stats_dt$P25[i], 0),
        format_number(stats_dt$Median[i], 0),
        format_number(stats_dt$P75[i], 0),
        format_number(stats_dt$Max[i], 0)
      )
    } else if (var_name %in% c("hourly_wage")) {
      row_values <- c(
        format_number(stats_dt$N[i], 0),
        format_number(stats_dt$Mean[i], 2),
        format_number(stats_dt$SD[i], 2),
        format_number(stats_dt$Min[i], 2),
        format_number(stats_dt$P25[i], 2),
        format_number(stats_dt$Median[i], 2),
        format_number(stats_dt$P75[i], 2),
        format_number(stats_dt$Max[i], 2)
      )
    } else if (var_name %in% c("wfh_0", "wfh_1", "wfh_2", "wfh_3", "wfh_4", "wfh_5", "salaried_exempt")) {
      row_values <- c(
        format_number(stats_dt$N[i], 0),
        format_number(stats_dt$Mean[i], 3),
        format_number(stats_dt$SD[i], 3),
        format_number(stats_dt$Min[i], 3),
        format_number(stats_dt$P25[i], 3),
        format_number(stats_dt$Median[i], 3),
        format_number(stats_dt$P75[i], 3),
        format_number(stats_dt$Max[i], 3)
      )
    } else {
      row_values <- c(
        format_number(stats_dt$N[i], 0),
        format_number(stats_dt$Mean[i], 1),
        format_number(stats_dt$SD[i], 1),
        format_number(stats_dt$Min[i], 1),
        format_number(stats_dt$P25[i], 1),
        format_number(stats_dt$Median[i], 1),
        format_number(stats_dt$P75[i], 1),
        format_number(stats_dt$Max[i], 1)
      )
    }
    
    latex_lines <- c(latex_lines,
                     paste(label, paste(row_values, collapse = " & "), "\\\\", sep = " & "))
  }
  
  # Close table
  latex_lines <- c(latex_lines,
                   "\\bottomrule",
                   "\\end{tabular}")
  
  return(paste(latex_lines, collapse = "\n"))
}

# Generate and save the LaTeX table
latex_table <- create_latex_table(summary_stats)
cat("\n================ SUMMARY STATISTICS TABLE (LaTeX) ================\n")
cat(latex_table)
cat("\n==================================================================\n")
cat(latex_table, file = file.path(tabs, "summary_stats_table.tex"))

# Summary by treatment status
if ("X" %in% names(dat)) {
  dat[, treatment := ifelse(X >= 0, "Above Threshold", "Below Threshold")]
  
  summary_by_treatment <- dat[, {
    calc_summary_stats(.SD, summary_vars)
  }, by = treatment]
  
  cat("\n================ SUMMARY BY TREATMENT STATUS ================\n")
  print(summary_by_treatment, digits = 3)
  cat("\n============================================================\n")
}




################################################################################
##  5. Compositional Data Analysis (Work from Home)
################################################################################

# Map every composition to the sphere
comp_cols <- paste0("wfh_", 0:5)

dat <- copy(temp)
dat[, wfh_vec := purrr::pmap(.SD, c), .SDcols = comp_cols]

eps <- 0.01
p   <- length(dat$wfh_vec[[1]])

dat[, wfh_vec_adj := lapply(wfh_vec, function(z) {
  z_eps <- (z + eps/p) / (1 + eps)
  y     <- sqrt(z_eps)
  y / sqrt(sum(y^2))
})]

dat[, wfh_sphere := lapply(wfh_vec, function(z) {
  y <- sqrt(z)
  y / sqrt(sum(y^2))
})]

# Run Fréchet test on compositional data
res_comp <- frechesTest(
  Y_obj    = dat$wfh_vec_adj,
  X_scalar = dat$X,
  c_val    = 0.0001,
  metric_space_type  = "sphere",
  h_frechet          = "CV",
  frechet_optns    = list(enforce_positive = TRUE),
  cv_K_folds         = 5,
  cv_n_bw_candidates = 10,
  verbose            = TRUE
)

cat(sprintf(
  "\n== Compositional Analysis ==\nN = %d  |  h = %.3f  |  T = %.3f  |  p = %.4f\n",
  nrow(dat), res_comp$h_mean_cv_selected, res_comp$Tn, res_comp$p_value))

# Helper: sphere to composition
to_comp <- function(y) {
  z <- y^2
  z / sum(z)
}

# Back-transform the Fréchet means
Z_left  <- to_comp(res_comp$l_hat_minus)
Z_right <- to_comp(res_comp$l_hat_plus)

# Create plot data
plot_dt <- data.table(
  side      = rep(c("Left of c", "Right of c"), each = length(comp_cols)),
  component = factor(rep(comp_cols, 2), levels = comp_cols),
  share     = c(Z_left, Z_right)
)

# Create better labels for the legend
wfh_labels <- c(
  "wfh_0" = "0 days (On-site)",
  "wfh_1" = "1 day",
  "wfh_2" = "2 days", 
  "wfh_3" = "3 days",
  "wfh_4" = "4 days",
  "wfh_5" = "5 days (Fully remote)"
)

# Reorder the factor levels
plot_dt[, component := factor(component, 
                              levels = comp_cols,
                              labels = wfh_labels[comp_cols])]

# Create the stacked bar chart
p_frechet <- ggplot(plot_dt, aes(x = side, y = share, fill = component)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     expand = c(0, 0),
                     limits = c(0, 1.02)) +
  scale_fill_manual(
    values = rev(project_colors$sequential[2:7]),
    name = "Work from Home\nDays per Week"
  ) +
  scale_x_discrete(labels = c("Left of c" = "Below Cutoff\n(Untreated)",
                              "Right of c" = "Above Cutoff\n(Treated)")) +
  labs(
    x = NULL, 
    y = "Share of Time"
  ) +
  theme_project() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(), 
    # the rectangle *around* the entire legend
    legend.background = element_rect(fill = "white", colour = NA),
    # the rectangles *behind* each key symbol
    legend.key        = element_rect(fill = "white", colour = NA)
  ) 

#--------------------------------------------
#####------- Figure 3 #####------- 
#--------------------------------------------
print(p_frechet)
save_plot(p_frechet, file.path("frechet_mean_composition"), width = 10, height = 5)

################################################################################
##  6. Compositional RD Plot with Binning
################################################################################

# Settings
bin_width <- 20000
cutoff    <- 0
comp_cols <- paste0("wfh_", 0:5)
x_limits     <- c(-100000, 100000)

# Create binned data
plot_dt <- copy(temp)
plot_dt[, bin_center := (floor(X / bin_width) + 0.5) * bin_width]

# Average composition within each bin
bin_means <- plot_dt[, lapply(.SD, mean, na.rm = TRUE),
                     by = bin_center,
                     .SDcols = comp_cols]

# Long format
long_dt_binned <- setDT(melt(
  bin_means,
  id.vars       = "bin_center",
  variable.name = "component",
  value.name    = "share"
))[, component := factor(component, levels = comp_cols)]

# Create binned plot without legend
p_binned <- ggplot(long_dt_binned, aes(x = bin_center, y = share, fill = component)) +
  geom_col(width = bin_width * 0.9, position = "fill") +
  geom_vline(xintercept = cutoff, linetype = "dashed", colour = "darkgrey", linewidth = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = rev(project_colors$sequential[2:7]),
    labels = rev(wfh_labels[comp_cols]),
    name = "Work from Home\nDays per Week"
  ) +
  labs(
    title = "A. Binned Composition",
    subtitle = sprintf("Bin width = %s", scales::dollar(bin_width)),
    x = NULL,
    y = "Share of Time"
  ) +
  theme_project() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = rel(1.1), face = "bold")
  ) + 
  coord_cartesian(xlim = x_limits)

################################################################################
##  7. Rolling Mean Compositional Plot
################################################################################

# Settings
window_width <- 20000
step_size    <- window_width/4
cutoff       <- 0
x_limits     <- c(-100000, 100000)

half_win <- window_width / 2
plot_dt  <- copy(temp)
setkeyv(plot_dt, "X")

# Helper: mean composition in one window
window_comp <- function(centre){
  idx <- plot_dt$X >= (centre - half_win) & plot_dt$X < (centre + half_win)
  if (!any(idx)) return(NULL)
  data.table(centre = centre,
             t(colMeans(plot_dt[idx, ..comp_cols], na.rm = TRUE)))
}

# Left-side windows
centres_left <- seq(
  from = floor((min(plot_dt$X) + half_win) / step_size) * step_size,
  to   = -half_win,
  by   = step_size)

left_dt <- rbindlist(
  Filter(Negate(is.null), lapply(centres_left, window_comp)),
  use.names = TRUE)
if (nrow(left_dt)) left_dt[, side := "Left"]

# Right-side windows
centres_right <- seq(
  from = half_win,
  to   = ceiling((max(plot_dt$X) - half_win) / step_size) * step_size,
  by   = step_size)

right_dt <- rbindlist(
  Filter(Negate(is.null), lapply(centres_right, window_comp)),
  use.names = TRUE)
if (nrow(right_dt)) right_dt[, side := "Right"]

# Combine and long format
roll_dt <- rbindlist(list(left_dt, right_dt), use.names = TRUE, fill = TRUE)

long_dt_rolling <- setDT(melt(
  roll_dt,
  id.vars       = c("side", "centre"),
  variable.name = "component",
  value.name    = "share"
))[, component := factor(component, levels = comp_cols)]

# Add connecting points at cutoff
left_last <- long_dt_rolling[side == "Left" & centre == max(centre[side == "Left"]), ]
left_last[, centre := cutoff]

right_first <- long_dt_rolling[side == "Right" & centre == min(centre[side == "Right"]), ]
right_first[, centre := cutoff]

long_dt_connected <- rbindlist(list(
  long_dt_rolling[side == "Left"],
  left_last,
  right_first,
  long_dt_rolling[side == "Right"]
))

# Create rolling mean plot with legend
p_rolling <- ggplot() +
  geom_area(data = long_dt_connected[side == "Left"],
            aes(x = centre, y = share, fill = component, group = component),
            position = "stack", alpha = 0.9) +
  geom_area(data = long_dt_connected[side == "Right"],
            aes(x = centre, y = share, fill = component, group = component),
            position = "stack", alpha = 0.9) +
  geom_vline(xintercept = cutoff, linetype = "dashed", colour = "darkgrey", linewidth = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = rev(project_colors$sequential[2:7]),
    labels = rev(wfh_labels[comp_cols]),
    name = "Work from Home\nDays per Week"
  ) +
  labs(
    title = "B. Rolling Mean Composition",
    subtitle = sprintf("Window = %s | Step = %s", dollar(window_width), dollar(step_size)),
    x = "Salary (centered at cutoff)",
    y = "Share of Time"
  ) +
  theme_project() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = rel(1.1), face = "bold"), 
    # the rectangle *around* the entire legend
    legend.background = element_rect(fill = "white", colour = NA),
    # the rectangles *behind* each key symbol
    legend.key        = element_rect(fill = "white", colour = NA)
  ) +
  coord_cartesian(xlim = x_limits)

# Combine plots using patchwork
library(patchwork)

# Extract legend from rolling plot
legend <- cowplot::get_legend(p_rolling)

# Remove legend from rolling plot for combined figure
p_rolling_no_legend <- p_rolling + theme(legend.position = "none")

# Combine plots
p_combined <- (p_binned / p_rolling_no_legend) + 
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = rel(1.3), face = "bold", hjust = 0),
      plot.subtitle = element_text(size = rel(1.05), hjust = 0),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# Add legend to the right
p_final <- cowplot::plot_grid(p_combined, legend, ncol = 2, rel_widths = c(3, 0.5))

#--------------------------------------------
#####------- Figure 2 #####------- 
#--------------------------------------------
# Print and save
print(p_final)
save_plot(p_final, file.path("wfh_composition_combined"), width = 10, height = 8)

################################################################################
##  8. Comparison of Estimation Methods
################################################################################

# Extract bandwidth
h_val <- res_comp$h_mean_cv_selected

# Subset to bandwidth
dat_band <- dat[abs(X) <= h_val]
dat_band[, side := ifelse(X < 0, "Left", "Right")]

# Triangular kernel weight
dat_band[, kern_w := pmax(0, 1 - abs(X)/h_val)]

# Flat band mean
flat_dt <- dat_band[, lapply(.SD, mean), by = side, .SDcols = comp_cols]
flat_dt[, Method := "Flat_Band"]

# Weighted plain mean
weighted_plain_dt <- dat_band[, {
  sum_w <- sum(kern_w)
  tmp <- lapply(.SD, function(z) sum(z * kern_w) / sum_w)
  as.list(tmp)
}, by = side, .SDcols = comp_cols]
weighted_plain_dt[, Method := "Weighted_Plain"]

# Spherical Fréchet mean
sph_list <- lapply(c("Left","Right"), function(s) {
  sub <- dat_band[side == s]
  Y_sph <- t(
    apply(sub[, ..comp_cols], 1L, function(z_vec) {
      y_raw  <- sqrt(z_vec)
      y_norm <- sqrt(sum(y_raw^2))
      y_raw / y_norm
    })
  )
  w_vec <- sub$kern_w
  sum_w  <- sum(w_vec)
  y_bar  <- colSums(Y_sph * w_vec) / sum_w
  y_bar_unit <- y_bar / sqrt(sum(y_bar^2))
  z_fr <- y_bar_unit^2
  as.list(setNames(z_fr, comp_cols))
})

sph_dt <- rbindlist(sph_list)
sph_dt[, side := c("Left","Right")]
sph_dt[, Method := "Sph_Frechet"]

# Combine all methods
compare_dt <- rbindlist(list(flat_dt, weighted_plain_dt, sph_dt), use.names=TRUE)
setcolorder(compare_dt, c("Method","side", comp_cols))

cat("\n--- Comparison of Local Mean Estimates (by side) ---\n")
print(compare_dt, digits = 4)

################################################################################
##  9. RD Diagnostics
################################################################################

dat_il <- copy(temp)

# McCrary density test
dens_out <- rddensity(dat_il$X, c = 0.000001, massPoints=TRUE)

cat("\n================  McCrary Density Test  ================\n")
print(summary(dens_out))


# Covariate continuity tests
dat_il$female <- as.numeric(dat_il$sex == 2)
covariates <- c("age_start","female","educ_level", "work_diff_state", "moved_state", 
                "work_diff_state_next", "moved_state_next",
                "wfh_any_now", "wfh_full_now", 
                "wfh_0", "wfh_1", "wfh_2", "wfh_3", "wfh_4", "wfh_5")

covariate_labels <- c("Age", "Female", "Education Level", 
                      "Work in Different State", "Moved State",
                      "Work in Different State Next", "Moved State Next",
                      "WFH Any", "WFH Full-Time",
                      "WFH 0 days/week", "WFH 1 day/week", 
                      "WFH 2 days/week", "WFH 3 days/week", 
                      "WFH 4 days/week", "WFH 5 days/week")

# Binary variables
dat_il[, wfh_full_next := as.numeric(wfh_full_next)]
dat_il[, wfh_full_now := as.numeric(wfh_full_now)]
dat_il[, moved_state := as.numeric(moved_state)]

# Balance tests
balance_tbl <- rbindlist(lapply(covariates, \(v) {
  tryCatch({
    w <- !is.na(dat_il[[v]])
    print(v)
    
    est <- rdrobust(y = dat_il[[v]][w],
                    x = dat_il$X[w],
                    c = 0.001, p = 1)
    
    test <- RDestimate(data = dat_il[w,], 
                       formula = as.formula(paste(v, "~ X")),
                       cutpoint = 0.0001, 
                       bw = est$bws['b', 'left'])
    
    # Create RD plot for each covariate
    p_rd <- rdplot(y = dat_il[[v]][w],
                   x = dat_il$X[w],
                   c = 0.0001, 
                   binselect = "esmv",
                   title = paste("RD Plot:", v),
                   x.label = "Running Variable",
                   y.label = v)
    
    data.table(
      covariate = v,
      tau_hat   = est$coef[1],
      se_tau    = est$se[1],
      p_value   = est$pv[[3]]
    )
  }, error = function(e) {
    warning(sprintf("Error for covariate %s: %s", v, e$message))
    data.table(
      covariate = v,
      tau_hat   = NA_real_,
      se_tau    = NA_real_,
      p_value   = NA_real_
    )
  })
}))

cat("\n================  RD Covariate Balance  ================\n")
print(balance_tbl, digits = 3)

# Print diagnostic info about p-values
cat("\n================  P-value Summary  ================\n")
cat("Number of coefficients with p < 0.01: ", sum(balance_tbl$p_value < 0.01, na.rm = TRUE), "\n")
cat("Number of coefficients with p < 0.05: ", sum(balance_tbl$p_value < 0.05, na.rm = TRUE), "\n")
cat("Number of coefficients with p < 0.10: ", sum(balance_tbl$p_value < 0.10, na.rm = TRUE), "\n")
cat("Range of p-values: ", min(balance_tbl$p_value, na.rm = TRUE), " to ", max(balance_tbl$p_value, na.rm = TRUE), "\n")

# Adjusted p-values
balance_tbl[, p_adj := p.adjust(p_value, method = "holm")]

cat("\n================  Results Summary  ================\n")
cat(sprintf("Fréchet test p-value: %.4f\n", res$p_value))
cat(sprintf("Number of observations: %d\n", nrow(dat)))
cat(sprintf("Bandwidth selected: %.0f\n", res$h_mean_cv_selected))
cat("\n")

################################################################################
##  12. Create LaTeX Table for RDD Results
################################################################################

# Function to format numbers for LaTeX
format_latex_number <- function(x, digits = 3) {
  if (is.na(x)) return("")
  formatted <- formatC(x, format = "f", digits = digits)
  # Remove trailing zeros after decimal point
  formatted <- gsub("(\\.\\d*?)0+$", "\\1", formatted)
  formatted <- gsub("\\.$", "", formatted)
  return(formatted)
}

# Create LaTeX table for balance tests
create_rdd_latex_table <- function(balance_tbl, covariate_labels = NULL) {
  # Add labels if provided before filtering
  if (!is.null(covariate_labels) && length(covariate_labels) == nrow(balance_tbl)) {
    balance_tbl$label <- covariate_labels
  } else {
    balance_tbl$label <- balance_tbl$covariate
  }
  
  # Filter out rows where coefficient is NA
  balance_tbl <- balance_tbl[!is.na(tau_hat)]
  
  # Start building LaTeX table
  latex_lines <- c(
    "\\begin{tabular}{lcccc}",
    "\\toprule",
    "Covariate & Estimate & Std. Error & p-value & Adj. p-value \\\\",
    "\\midrule"
  )
  
  # Add data rows
  for (i in 1:nrow(balance_tbl)) {
    # Format estimate with stars for significance
    est_str <- format_latex_number(balance_tbl$tau_hat[i], 3)
    
    # Add significance stars based on p-value
    if (!is.na(balance_tbl$p_value[i])) {
      if (balance_tbl$p_value[i] < 0.01) {
        est_str <- paste0(est_str, "***")
      } else if (balance_tbl$p_value[i] < 0.05) {
        est_str <- paste0(est_str, "**")
      } else if (balance_tbl$p_value[i] < 0.10) {
        est_str <- paste0(est_str, "*")
      }
    }
    
    # Format other values
    se_str <- paste0("(", format_latex_number(balance_tbl$se_tau[i], 3), ")")
    p_str <- format_latex_number(balance_tbl$p_value[i], 3)
    p_adj_str <- format_latex_number(balance_tbl$p_adj[i], 3)
    
    # Create row
    row_str <- paste(balance_tbl$label[i], est_str, se_str, p_str, p_adj_str, sep = " & ")
    latex_lines <- c(latex_lines, paste0(row_str, " \\\\"))
  }
  
  # Add bottom of table
  latex_lines <- c(
    latex_lines,
    "\\midrule",
    paste0("N & \\multicolumn{4}{c}{", format(nrow(dat_il), big.mark = ","), "} \\\\"),
    "\\bottomrule",
    "\\end{tabular}"
  )
  
  return(paste(latex_lines, collapse = "\n"))
}

# Create the LaTeX table
latex_rdd_table <- create_rdd_latex_table(balance_tbl, covariate_labels)

# Print to console
cat("\n================ RDD BALANCE TESTS TABLE (LaTeX) ================\n")
cat(latex_rdd_table)
cat("\n=================================================================\n")

# Save to file
writeLines(latex_rdd_table, file.path(tabs, "rdd_balance_tests.tex"))

# Also create a more detailed regression table with additional statistics
create_detailed_rdd_latex <- function(balance_tbl, covariate_labels = NULL) {
  # Add labels if provided before filtering
  if (!is.null(covariate_labels) && length(covariate_labels) == nrow(balance_tbl)) {
    balance_tbl$label <- covariate_labels
  } else {
    balance_tbl$label <- balance_tbl$covariate
  }
  
  # Filter out rows where coefficient is NA
  balance_tbl <- balance_tbl[!is.na(tau_hat)]
  
  # Start building LaTeX table
  latex_lines <- c(
    "\\begin{tabular}{lc}",
    "\\toprule",
    "Covariate & Coefficient \\\\",
    " & (Std. Error) \\\\",
    "\\midrule"
  )
  
  # Add data rows
  for (i in 1:nrow(balance_tbl)) {
    # Format estimate with stars
    est_str <- format_latex_number(balance_tbl$tau_hat[i], 3)
    
    # Add significance stars
    if (!is.na(balance_tbl$p_value[i])) {
      if (balance_tbl$p_value[i] < 0.01) {
        est_str <- paste0(est_str, "***")
      } else if (balance_tbl$p_value[i] < 0.05) {
        est_str <- paste0(est_str, "**")
      } else if (balance_tbl$p_value[i] < 0.10) {
        est_str <- paste0(est_str, "*")
      }
    }
    
    # Format standard error
    se_str <- paste0("(", format_latex_number(balance_tbl$se_tau[i], 3), ")")
    
    # Create main row and sub-row
    latex_lines <- c(
      latex_lines,
      paste0(balance_tbl$label[i], " & ", est_str, " \\\\"),
      paste0(" & ", se_str, " \\\\")
    )
    
    # Add spacing between variables (except after last one)
    if (i < nrow(balance_tbl)) {
      latex_lines <- c(latex_lines, "[0.5ex]")
    }
  }
  
  # Add bottom of table
  latex_lines <- c(
    latex_lines,
    "\\midrule",
    paste0("Observations & ", format(nrow(dat_il), big.mark = ","), " \\\\"),
    "\\bottomrule",
    "\\end{tabular}"
  )
  
  return(paste(latex_lines, collapse = "\n"))
}

# Create detailed table
latex_detailed_table <- create_detailed_rdd_latex(balance_tbl, covariate_labels)

# Save detailed version
writeLines(latex_detailed_table, file.path(tabs, "rdd_balance_detailed.tex"))

cat("LaTeX tables saved to:")
cat("\n  - ", file.path(tabs, "rdd_balance_tests.tex"))
cat("\n  - ", file.path(tabs, "rdd_balance_detailed.tex"))
