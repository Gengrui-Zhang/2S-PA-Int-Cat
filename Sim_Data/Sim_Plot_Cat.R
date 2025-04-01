library(here)
library(dplyr)
library(tidyr)
library(forcats)
library(readr)
library(ggplot2)
library(rlang)

# Read Data
data_path <- here("Sim_Data", "categorical_03252025.rds")
df <- readRDS(data_path)

# Helper Function
process_data <- function(data) {
  results <- data %>%
    gather("var", "val", raw_bias.upi.xm_est_std:rmse.lmsfs.xm_est_std) %>%
    dplyr::select(-c(REPLICATIONS:WARNINGS)) %>%
    separate(col = var, into = c("stats", "method"), sep = "\\.") %>%
    spread(stats, val) %>%
    mutate(N_lab = as_factor(paste0("italic(N) == ", N)),
           rel_lab = as_factor(paste0("rho == ", rel)),
           gammaint_lab = as_factor(paste0("\\beta_{xm} == ", gamma_xm)),
           skew_lab = as_factor(paste0("Skewness == ", skewness)))
  
  return(results)
}

# Plot function
plot_results <- function(data,
                         x_var, y_var, 
                         x_label, y_label, 
                         y_limits = NULL, 
                         color_var = "method", 
                         shape_var = "method", 
                         point_size = 3, 
                         add_lines = FALSE, 
                         remove_baseline = FALSE,
                         method_labels = NULL,
                         method_colors = NULL,
                         method_shapes = NULL,
                         method_linetypes = NULL) {
  
  # Extract unique methods if not provided
  unique_methods <- sort(unique(data[[color_var]]))
  
  # Set default labels, colors, shapes, linetypes if not specified
  if (is.null(method_labels)) method_labels <- setNames(unique_methods, unique_methods)
  if (is.null(method_colors)) method_colors <- setNames(RColorBrewer::brewer.pal(length(unique_methods), "Set2"), unique_methods)
  if (is.null(method_shapes)) method_shapes <- setNames(seq(15, 15 + length(unique_methods) - 1), unique_methods)
  if (is.null(method_linetypes)) method_linetypes <- setNames(rep("solid", length(unique_methods)), unique_methods)
  
  plot <- data %>%
    ggplot(aes(x = factor(!!sym(x_var)), y = !!sym(y_var), color = !!sym(color_var), shape = !!sym(shape_var))) +
    geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0), size = point_size, alpha = 0.8) +
    facet_grid(skew_lab ~ rel_lab, labeller = label_parsed) +
    scale_shape_manual(name = "Method", values = method_shapes, labels = method_labels) +
    scale_color_manual(name = "Method", values = method_colors, labels = method_labels) +
    labs(x = x_label, y = y_label, color = "Method", shape = "Method") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      strip.text = element_text(size = 12, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      strip.background = element_rect(color = "black", size = 1),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray80"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', colour = "gray80"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  
  if (add_lines) {
    plot <- plot + 
      geom_line(aes(group = !!sym(color_var), linetype = !!sym(shape_var)), size = 0.7, alpha = 0.6) +
      scale_linetype_manual(name = "Method", values = method_linetypes, labels = method_labels)
  }
  
  if (!is.null(y_limits)) {
    plot <- plot + ylim(y_limits)
  }
  
  return(plot)
}

# Split the data
type1data <- df %>% filter(gamma_xm == 0)
powerdata <- df %>% filter(gamma_xm == 0.3)

type1_pd <- process_data(type1data)
power_pd <- process_data(powerdata)

# Write the processed data
write_csv(type1_pd, "Sim_Data/type1_03252025.csv")
write_csv(power_pd, "Sim_Data/power_03252025.csv")

# Plot the results
type1_plot <- read.csv("Sim_Data/type1_03252025.csv")
# type1_plot <- type1_plot %>% filter(method != "reg")
power_plot <- read.csv("Sim_Data/power_03252025.csv")
# power_plot <- power_plot %>% filter(method != "reg")

# Standard Bias
methods <- c("upi", "tspa", "lms", "lmsfs")
method_labels = c(upi = "All-Pair UPI", tspa = "2S-PA-Int", 
                  lms = "LMS", lmsfs = "LMS_FS")
method_colors = c(upi = "#009E73", tspa = "#CC79A7", 
                  lms = "#E69F00", lmsfs = "pink") 
method_shapes = c(upi = 15, tspa = 18, 
                  lms = 16, lmsfs = 17)
method_linetypes = c(upi = "twodash", tspa = "dashed",
                     lms = "solid", lmsfs = "dotted")

# Raw Bias
raw_bias_plot <- plot_results(
  data = power_plot,
  x_var = "N",
  y_var = "raw_bias",
  x_label = "Number of Items",
  y_label = "Raw Bias",
  add_lines = TRUE,
  method_labels = method_labels,
  method_colors = method_colors,
  method_shapes = method_shapes,
  method_linetypes = method_linetypes
)

raw_bias_plot <- raw_bias_plot +
  geom_hline(yintercept = c(-0.4, 0.4), linetype = "dashed", color = "red")

# Standardized Bias
sd_bias_plot <- plot_results(
  data = type1_plot,
  x_var = "N",
  y_var = "std_bias",
  x_label = "Number of Items",
  y_label = "Standardized Bias",
  add_lines = TRUE,
  method_labels = method_labels,
  method_colors = method_colors,
  method_shapes = method_shapes,
  method_linetypes = method_linetypes
)

sd_bias_plot <- sd_bias_plot +
  geom_hline(yintercept = c(-0.4, 0.4), linetype = "dashed", color = "red")

# Median-Mad Relative SE Bias
stdMed_rse_bias_plot <- plot_results(
  data = power_plot,
  x_var = "N",
  y_var = "stdMed_rse_bias",
  x_label = "Number of Items",
  y_label = "Median-Mad Relative SE Bias",
  add_lines = TRUE,
  method_labels = method_labels,
  method_colors = method_colors,
  method_shapes = method_shapes,
  method_linetypes = method_linetypes
)

stdMed_rse_bias_plot <- stdMed_rse_bias_plot +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "red")

# Median-Mad Relative SE Bias
raw_bias_plot <- plot_results(
  data = power_plot,
  x_var = "N",
  y_var = "raw_rse_bias",
  x_label = "Number of Items",
  y_label = "Raw SE Bias",
  add_lines = TRUE,
  method_labels = method_labels,
  method_colors = method_colors,
  method_shapes = method_shapes,
  method_linetypes = method_linetypes
)

raw_bias_plot <- raw_bias_plot +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "red")

# Coverage rate
cov_plot <- plot_results(
  data = power_plot,
  x_var = "N",
  y_var = "coverage",
  x_label = "Number of Items",
  y_label = "Coverage",
  add_lines = TRUE,
  method_labels = method_labels,
  method_colors = method_colors,
  method_shapes = method_shapes,
  method_linetypes = method_linetypes
)

cov_plot <- cov_plot +
  geom_hline(yintercept = c(0.91), linetype = "dashed", color = "red")

# RMSE
rmse_plot <- plot_results(
  data = power_plot,
  x_var = "N",
  y_var = "rmse",
  x_label = "Number of Items",
  y_label = "RMSE",
  add_lines = TRUE,
  method_labels = method_labels,
  method_colors = method_colors,
  method_shapes = method_shapes,
  method_linetypes = method_linetypes
)

# Type I error rate
type1_plot <- plot_results(
  data = type1_plot,
  x_var = "N",
  y_var = "type1_std",
  x_label = "Number of Items",
  y_label = "Type I Error",
  add_lines = TRUE,
  method_labels = method_labels,
  method_colors = method_colors,
  method_shapes = method_shapes,
  method_linetypes = method_linetypes
)

# Statistical Power
power_plot <- plot_results(
  data = power_plot,
  x_var = "N",
  y_var = "power_std",
  x_label = "Number of Items",
  y_label = "Power",
  add_lines = TRUE,
  method_labels = method_labels,
  method_colors = method_colors,
  method_shapes = method_shapes,
  method_linetypes = method_linetypes
)
power_plot <- power_plot +
  geom_hline(yintercept = c(0.8), linetype = "dashed", color = "red")
