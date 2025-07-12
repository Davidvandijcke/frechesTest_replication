# ================================================================================
# REPLOT FRÉCHET ANOVA POWER CURVES FROM SAVED RESULTS
# 
# This script allows you to recreate the power curve plots without rerunning
# the simulations. You can also modify plot aesthetics as needed.
# ================================================================================

# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)

# Load the saved results
cat("Loading saved simulation results...\n")
all_power_results_df <- readRDS("../figs/raw_results/frechetANOVA_power_results.rds")
simulation_params <- readRDS("../figs/raw_results/frechetANOVA_simulation_params.rds")

# Restore key parameters
list2env(simulation_params, envir = .GlobalEnv)

cat("Loaded results for N =", SAMPLE_SIZE_POWER_CURVE, "with", 
    nrow(all_power_results_df), "parameter configurations\n\n")

# --- RECREATE PLOT WITH CUSTOM MODIFICATIONS ---

# You can modify these aesthetics as needed
test_colors <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1:2]
# Or use custom colors:
# test_colors <- c("Fréchet Jump Test" = "#E41A1C", "Fréchet ANOVA" = "#377EB8")

test_shapes <- c(16, 17)
test_linetypes <- c("solid", "longdash")
test_labels <- c("Fréchet Jump Test", "Fréchet ANOVA")

names(test_colors) <- test_labels
names(test_shapes) <- test_labels
names(test_linetypes) <- test_labels

# Initialize plot list
plot_list_final <- list()
subplot_idx_counter <- 0

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
    subplot_title_text_full <- paste0("(", LETTERS[subplot_idx_counter], ") ", 
                                      tools::toTitleCase(metric_plot_loop), 
                                      " (", tools::toTitleCase(dgp_plot_loop), " DGP)")
    
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
    
    # Create subplot with customizable text sizes
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
      theme_classic(base_size = 14) +  # Adjust base size as needed
      theme(
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = rel(1.1), face = "plain"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"), 
        axis.title.x = element_text(size = rel(1.0), margin = margin(t = 5)),
        axis.title.y = element_text(size = rel(1.0), margin = margin(r = 5)),
        axis.text = element_text(size = rel(0.9))
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

# Combine plots in 2x3 layout
if (length(plot_list_final) == 6) {
  combined_figure_patch <- (plot_list_final[[1]] | plot_list_final[[2]] | plot_list_final[[3]]) / 
                           (plot_list_final[[4]] | plot_list_final[[5]] | plot_list_final[[6]])
  
  # Add legend and caption with clean styling
  final_figure_with_legend <- combined_figure_patch + 
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = paste0("N = ", SAMPLE_SIZE_POWER_CURVE, 
                       ", Test Level α = ", ALPHA_LEVEL,
                       ". Results based on ", N_SIMULATIONS_POWER_CURVE, 
                       " Monte Carlo simulations per point.")
    ) &
    theme(
      legend.position = "bottom", 
      legend.box.margin = margin(t = 15, b = 5), 
      legend.title = element_text(face = "bold", size = rel(1.2)),
      legend.text = element_text(size = rel(1.1)),
      plot.caption = element_text(hjust = 0.5, size = rel(1.0),
                                  margin = margin(t = 15, b = 5)),
      # Clean legend (no grey box)
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA)
    )
  
  # Display plot
  print(final_figure_with_legend)
  
  # Save plot with custom filename if desired
  output_filename <- paste0("../figs/power_curves_stacked_2x3_N", 
                            SAMPLE_SIZE_POWER_CURVE, "_replot.png")
  
  ggsave(output_filename, 
         plot = final_figure_with_legend, 
         width = 12, height = 8, dpi = 300,
         bg = "white")
  
  cat(paste0("\nPlot saved to: ", output_filename, "\n"))
  
  # Optional: Save as PDF for publication
  # ggsave(gsub(".png", ".pdf", output_filename), 
  #        plot = final_figure_with_legend, 
  #        width = 12, height = 8,
  #        device = cairo_pdf)
  
} else {
  cat("Error: Incorrect number of plots generated.\n")
}

# --- ADDITIONAL PLOT OPTIONS ---

# Example: Create individual plots for each metric space
# for (metric in unique(all_power_results_df$metric_space)) {
#   subset_data <- all_power_results_df %>% filter(metric_space == metric)
#   # Create plot for this metric...
# }

# Example: Change to 3x2 layout (original)
# if (length(plot_list_final) == 6) {
#   combined_figure_3x2 <- (plot_list_final[[1]] | plot_list_final[[4]]) / 
#                          (plot_list_final[[2]] | plot_list_final[[5]]) /
#                          (plot_list_final[[3]] | plot_list_final[[6]])
#   # Continue with legend...
# }