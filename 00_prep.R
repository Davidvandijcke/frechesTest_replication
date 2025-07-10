
#****************************************************************************************************************************************************

# LOAD LIBRARIES

#****************************************************************************************************************************************************

# load packages
packages_not_load <- c()
packages_load <- c(
  # Data manipulation
  "data.table",
  
  # Statistical analysis
  "rddisttest",
  "frechesTest",
  "rdd",
  "rddensity",
  "rdrobust",
  
  # Visualization
  "ggplot2",
  "scales",
  "viridis",
  "RColorBrewer",
  "reshape2",
  
  # Utilities
  "purrr"
)


{
sink("/dev/null") # load packages but suppress output
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = packages_load, character.only = TRUE)


sink()

}

source("utils.R")

#****************************************************************************************************************************************************

# SET AESTHETIC THEME AND COLOR SCHEME

#****************************************************************************************************************************************************

# Define custom color scheme - modern vibrant palette
# Option 1: Vibrant Modern (current)
project_colors <- list(
  # Main palette - vibrant modern colors
  main = c("#FF006E", "#8338EC", "#3A86FF", "#06FFB4", "#FB5607", "#FFBE0B"),
  
  # Diverging palette for heatmaps (pink to teal through white)
  diverging = c("#FF006E", "#F72585", "#E85D75", "#FDDBC7", "#C8E7E1", "#06FFB4", "#00D9A3", "#00BF93"),
  
  # Sequential palette (sunset gradient)
  sequential = c("#FFF3E0", "#FFE0B2", "#FFCC80", "#FFB74D", "#FFA726", "#FF9800", "#FB8C00", "#F57C00", "#EF6C00"),
  
  # Categorical palette (modern vibrant mix)
  categorical = c("#FF006E", "#8338EC", "#3A86FF", "#06FFB4", "#FB5607", "#FFBE0B", "#00BF63", "#C77DFF", "#7209B7", "#560BAD")
)

# Option 2: Miami Vice / Synthwave
project_colors <- list(
  main = c("#FF006E", "#00D9FF", "#FF5E00", "#00FF88", "#6E00FF", "#FFD600"),
  diverging = c("#FF006E", "#FF4081", "#FF80AB", "#FFE4EC", "#E0F7FA", "#80DEEA", "#26C6DA", "#00ACC1"),
  sequential = c("#FCE4EC", "#F8BBD0", "#F48FB1", "#F06292", "#EC407A", "#E91E63", "#D81B60", "#C2185B", "#AD1457"),
  categorical = c("#FF006E", "#00D9FF", "#FF5E00", "#00FF88", "#6E00FF", "#FFD600", "#FF1744", "#00E5FF", "#FF6E40", "#69F0AE")
)

# Option 3: Pastel Dream
project_colors <- list(
  main = c("#FFB3BA", "#BAFFC9", "#BAE1FF", "#FFFFBA", "#FFDFBA", "#E0BBE4"),
  diverging = c("#FFB3BA", "#FFC4C9", "#FFD5D8", "#FFE6E7", "#E7F3FF", "#D1E7FF", "#BAE1FF", "#A3D5FF"),
  sequential = c("#FFF5F5", "#FFE0E6", "#FFCCD5", "#FFB3BA", "#FF99A1", "#FF8089", "#FF6670", "#FF4D58", "#FF3340"),
  categorical = c("#FFB3BA", "#BAFFC9", "#BAE1FF", "#FFFFBA", "#FFDFBA", "#E0BBE4", "#D4A5A5", "#A8E6CF", "#C7CEEA", "#FFDAC1")
)

# Option 4: Blue-Gray Gradient
project_colors <- list(
  main = c("#F0F2F8FF", "#D2DCF0FF", "#A0B7D8FF", "#7696BEFF", "#606A81FF", "#384157FF"),
  diverging = c("#384157FF", "#606A81FF", "#7696BEFF", "#A0B7D8FF", "#D2DCF0FF", "#F0F2F8FF", "#F0F2F8FF", "#FFFFFF"),
  sequential = c("#F0F2F8FF", "#D2DCF0FF", "#A0B7D8FF", "#7696BEFF", "#606A81FF", "#384157FF"),
  categorical = c("#F0F2F8FF", "#D2DCF0FF", "#A0B7D8FF", "#7696BEFF", "#606A81FF", "#384157FF")
)

# Define consistent theme for all plots - modern and clean
theme_project <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Text elements - darker and more contrast
      plot.title = element_text(size = rel(1.3), face = "bold", color = "#1a1a1a", hjust = 0, margin = margin(b = 12)),
      plot.subtitle = element_text(size = rel(1.05), color = "#404040", hjust = 0, margin = margin(b = 15)),
      plot.caption = element_text(size = rel(0.85), color = "#606060", hjust = 1, margin = margin(t = 15)),
      
      # Axis elements - cleaner look
      axis.title = element_text(size = rel(1.05), face = "bold", color = "#1a1a1a"),
      axis.title.x = element_text(margin = margin(t = 12)),
      axis.title.y = element_text(margin = margin(r = 12)),
      axis.text = element_text(size = rel(0.95), color = "#2a2a2a"),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "#cccccc", linewidth = 0.4),
      
      # Legend elements - more prominent
      legend.title = element_text(size = rel(1.05), face = "bold", color = "#1a1a1a"),
      legend.text = element_text(size = rel(0.95), color = "#2a2a2a"),
      legend.key.size = unit(0.9, "cm"),
      legend.background = element_rect(fill = "#fafafa", color = NA),
      legend.box.background = element_rect(fill = "#fafafa", color = NA),
      
      # Panel elements - lighter and cleaner
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      
      # Plot margins and background
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Strip text for facets - modern style
      strip.background = element_rect(fill = "#f5f5f5", color = NA),
      strip.text = element_text(size = rel(1), face = "bold", color = "#1a1a1a", margin = margin(8, 8, 8, 8))
    )
}

# Set the custom theme as default
theme_set(theme_project())

# Update default ggplot2 colors
options(
  ggplot2.discrete.colour = project_colors$categorical,
  ggplot2.discrete.fill = project_colors$categorical
)

# Utility functions
format_number <- function(x, digits = 2, big.mark = ",") {
  formatC(x, format = "f", digits = digits, big.mark = big.mark)
}

save_plot <- function(plot, filename, path = figs, width = 8, height = 6, dpi = 300) {
  base_path <- file.path(path, filename)
  
  # Save as PDF
  ggsave(
    filename = paste0(base_path, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = "pdf"
  )
  
  # Save as PNG
  ggsave(
    filename = paste0(base_path, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    device = "png"
  )
}
