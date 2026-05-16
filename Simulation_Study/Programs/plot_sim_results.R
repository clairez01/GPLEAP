# ------------------------------------------------------------------------------
# Program name: plot_sim_results.R
# Created by..: Claire Zhu
# Date created: 12 FEB 2026
# Description.: This program contains the functions for plotting the metrics
#               (bias, MSE, 95% CI coverage probability, and 95% CI width) in a 
#               simulation study
# ------------------------------------------------------------------------------

# load libraries
library(ggplot2)
library(patchwork)
library(ggpubr)
library(haven)
library(dplyr)

# 1. Data Cleaning and Transformation
# -----------------------------------------------------------------------------------
setwd("/nas/longleaf/home/clairez1/Paper2/Results")
df_gpleap  <- readRDS("compiled_results_gpleap.rds") 
df_gp      <- readRDS("compiled_results_gp.rds")
df_comp    <- readRDS("compiled_results_leap_pp.rds")
#df_gpleap  <- readRDS("compiled_results_gpleap_sens.rds") 
#df_gp      <- readRDS("compiled_results_gp_sens.rds")
#df_comp    <- readRDS("compiled_results_leap_pp_sens.rds")

# Run this for sensitivity analysis
#n0.filter <- 50
#df_gpleap <- subset(df_gpleap, n0 == n0.filter)
#df_gp     <- subset(df_gp, n0 == n0.filter)
#df_comp   <- subset(df_comp, n0 == n0.filter)


# --------------
# Pre-processing
#---------------
# Ensure numeric types
cols_to_fix <- c("bias", "MSE", "cov_prob", "width", "q", "pr")
df_gpleap[cols_to_fix] <- lapply(df_gpleap[cols_to_fix], function(x) as.numeric(as.character(x)))

# Remove the gamma parameter for the TE plots
df_plot <- df_gpleap %>% filter(method != "GPLEAP_gamma") %>% filter(method != "GPLEAP_Pooled")

# --- DATA DUPLICATION TRICK START ---
# Isolate the scenario where pr = 1 (where q is usually 0)
pr_one <- df_plot %>% filter(pr == 1)

# Get all unique q values from the rest of the simulation
all_q <- unique(df_plot$q[df_plot$q != 0])

# Duplicate the pr=1 results for every non-zero q value
# This creates a mathematical endpoint for every line in the plot
extra_points <- do.call(rbind, lapply(all_q, function(val) {
  temp <- pr_one
  temp$q <- val
  return(temp)
}))

# Combine original data with the duplicated points
df_plot <- rbind(df_plot, extra_points)
# --- DATA DUPLICATION TRICK END ---

# Create 'Proportion Non-exchangeable' 
df_plot$pr.u <- 1 - df_plot$pr

# Separate by data type (dt) - these will now contain the duplicated endpoints
df.linear <- df_plot %>% filter(dt == 1)
df.quad   <- df_plot %>% filter(dt == 2)
df.emax   <- df_plot %>% filter(dt == 3)

# Ensure numeric types for df_comp
cols_to_fix <- c("bias", "MSE", "cov_prob", "width", "q", "pr")
df_comp[cols_to_fix] <- lapply(df_comp[cols_to_fix], function(x) as.numeric(as.character(x)))

# --- DATA DUPLICATION TRICK FOR df_comp START ---
# Isolate the scenario where pr = 1
comp_pr_one <- df_comp %>% filter(pr == 1)

# Get all unique q values from the rest of the simulation
comp_all_q <- unique(df_comp$q[df_comp$q != 0])

# Duplicate the pr=1 results for every non-zero q value
comp_extra_points <- do.call(rbind, lapply(comp_all_q, function(val) {
  temp <- comp_pr_one
  temp$q <- val
  return(temp)
}))

# Combine original data with the duplicated points
df_comp <- rbind(df_comp, comp_extra_points)
# --- DATA DUPLICATION TRICK FOR df_comp END ---

# Create 'Proportion Non-exchangeable' 
df_comp$pr.u <- 1 - df_comp$pr

# Separate by data type (dt) - these will now contain the duplicated endpoints
df_comp.linear <- df_comp %>% filter(dt == 1)
df_comp.quad   <- df_comp %>% filter(dt == 2)
df_comp.emax   <- df_comp %>% filter(dt == 3)


# 2. Theme and Colors
# -----------------------------------------------------------------------------------
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
tableau_colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2") 

# 3. Create Plotting Function
# -----------------------------------------------------------------------------------
create_panels_lowdose <- function(data, title_suffix, bias, mse, cov, wid) {
  
  # 1. Define descriptive labels for the legend
  method_labels <- c(
    "LEAP" = "LEAP-linear", 
    "GPLEAP_Low" = "GPLEAP",
    "pp_0.5" = "PP-linear"
  )
  
  # 2. Absolute Bias Plot
  p_bias <- ggplot(data, aes(x = pr, y = abs(bias), color = method, linetype = factor(q))) +
    geom_point() + geom_line() +
    geom_hline(yintercept = bias[1], linewidth = .7, linetype = "solid", color = "black") +
    geom_hline(yintercept = bias[2], linewidth = .7, linetype = "solid", color = "#A6611D") +
    labs(title = paste("Absolute Bias:", title_suffix), 
         x = "Proportion Exchangeable", 
         y = "Absolute Bias") +
    scale_color_manual(values = c("#F28E2B", "#E15759", "#76B7B2"), labels = method_labels, name = "Method") +
    scale_linetype_discrete(name = expression(paste("Historical Bias (", delta, ")")))
  
  # 3. Relative MSE Plot
  p_mse_rel <- ggplot(data, aes(x = pr, y = MSE / mse[1], color = method, linetype = factor(q))) +
    geom_point() + geom_line() +
    geom_hline(yintercept = 1, linewidth = .8, linetype = "solid", color = "#2B435C") +
    geom_hline(yintercept = mse[2] / mse[1], linewidth = .8, linetype = "solid", color = "#A6611D") +
    labs(title = paste("Relative MSE:", title_suffix), 
         x = "Proportion Exchangeable", 
         y = "Relative MSE (vs ANOVA)") +
    scale_color_manual(values = c("#F28E2B", "#E15759", "#76B7B2"), labels = method_labels, name = "Method") +
    scale_linetype_discrete(name = expression(paste("Historical Bias (", delta, ")"))) +
    coord_cartesian(ylim = c(NA, 1.5)) # Limits top to 1.5x ANOVA, lets bottom auto-scale
  
  # 3. MSE Plot
  p_mse <- ggplot(data, aes(x = pr, y = MSE, color = method, linetype = factor(q))) +
    geom_point() + geom_line() +
    geom_hline(yintercept = mse[1], linewidth = .7, linetype = "solid", color = "black") +
    geom_hline(yintercept = mse[2], linewidth = .7, linetype = "solid", color = "#A6611D") +
    labs(title = paste("MSE:", title_suffix), 
         x = "Proportion Exchangeable", 
         y = "MSE") +
    scale_color_manual(values = c("#F28E2B", "#E15759", "#76B7B2"), labels = method_labels, name = "Method") +
    scale_linetype_discrete(name = expression(paste("Historical Bias (", delta, ")")))
  
  # 4. Coverage Plot
  p_cov <- ggplot(data, aes(x = pr, y = cov_prob, color = method, linetype = factor(q))) +
    geom_point() + geom_line() +
    geom_hline(yintercept = cov[1], linewidth = .7, linetype = "solid", color = "black") +
    geom_hline(yintercept = cov[2], linewidth = .7, linetype = "solid", color = "#A6611D") +
    labs(title = paste("Coverage:", title_suffix), 
         x = "Proportion Exchangeable", 
         y = "Coverage Probability") +
    scale_color_manual(values = c("#F28E2B", "#E15759", "#76B7B2"), labels = method_labels, name = "Method") +
    scale_linetype_discrete(name = expression(paste("Historical Bias (", delta, ")"))) +
    coord_cartesian(ylim = c(0.75, 1))
  
  # 5. Relative CrI Width Plot
  p_wid_rel <- ggplot(data, aes(x = pr, y = width / wid[1], color = method, linetype = factor(q))) +
    geom_point() + geom_line() +
    geom_hline(yintercept = 1, linewidth = .8, linetype = "solid", color = "#2B435C") +
    geom_hline(yintercept = wid[2] / wid[1], linewidth = .8, linetype = "solid", color = "#A6611D") +
    labs(title = paste("Relative CrI Width:", title_suffix), 
         x = "Proportion Exchangeable", 
         y = "Relative CrI Width (vs ANOVA)") +
    scale_color_manual(values = c("#F28E2B", "#E15759", "#76B7B2"), labels = method_labels, name = "Method") +
    scale_linetype_discrete(name = expression(paste("Historical Bias (", delta, ")"))) 
    #coord_cartesian(ylim = c(NA, 1.5))
  
  # 5. Width Plot
  p_wid <- ggplot(data, aes(x = pr, y = width, color = method, linetype = factor(q))) +
    geom_point() + geom_line() +
    geom_hline(yintercept = wid[1], linewidth = .7, linetype = "solid", color = "black") +
    geom_hline(yintercept = wid[2], linewidth = .7, linetype = "solid", color = "#A6611D") +
    labs(title = paste("CI Width:", title_suffix), 
         x = "Proportion Exchangeable", 
         y = "95% Credible Interval Width") +
    scale_color_manual(values = c("#F28E2B", "#E15759", "#76B7B2"), labels = method_labels, name = "Method") +
    scale_linetype_discrete(name = expression(paste("Historical Bias (", delta, ")")))
  
  # Arrange in 2x2
  ggarrange(p_bias, p_mse_rel, p_cov, p_wid_rel, ncol=2, nrow=2, 
            common.legend = TRUE, legend="bottom")
}

# 4. Generate the Final Plots
# -----------------------------------------------------------------------------------

# Plots for main paper

# ANOVA Bias
bias_lowdose <- c(2.2163057, 0.6123392, 1.0949452)
bias_highdose <- c(1.8225415, -2.4316934, 0.5147520)

# ANOVA MSE
mse_lowdose <- c(32.24274, 26.86635, 27.44296)
mse_highdose <- c(36.98988, 93.44788, 29.05168)

# ANOVA 95% CI Coverage
cov_lowdose <- c(0.9344, 0.9493, 0.9440)
cov_highdose <- c(0.9319, 0.8543, 0.9494)

# ANOVA Mean CI Width
width_lowdose <- c(21.04972, 26.09483, 20.14880)
width_highdose <- c(22.13176, 27.33523, 21.19022)

df_gpleap.linear.low <- df.linear %>% filter(method == "GPLEAP_Low")
df_gpleap.quad.low   <- df.quad %>% filter(method == "GPLEAP_Low")
df_gpleap.emax.low   <- df.emax %>% filter(method == "GPLEAP_Low")

df_comp.linear.low <- df_comp.linear %>% filter(method == "LEAP" | method == "pp_0.5")
df_comp.quad.low   <- df_comp.quad %>% filter(method == "LEAP" | method == "pp_0.5")
df_comp.emax.low   <- df_comp.emax %>% filter(method == "LEAP" | method == "pp_0.5")

df.linear.low <- rbind(df_gpleap.linear.low, df_comp.linear.low)
df.quad.low <- rbind(df_gpleap.quad.low, df_comp.quad.low)
df.emax.low <- rbind(df_gpleap.emax.low, df_comp.emax.low)

create_panels_lowdose(df.linear.low, "Linear", bias = abs(c(2.2163057, -0.57)), mse = c(32.24274, 17), cov = c(0.9344, 0.957), wid = c(21.04972, 16.7))
create_panels_lowdose(df.quad.low, "Quadratic", bias = abs(c(0.6123392, -2.21)), mse = c(27.86635, 28.1), cov = c(0.9493, 0.916), wid = c(26.09483, 18.6))
create_panels_lowdose(df.emax.low, "Emax", bias = abs(c(1.0949452, -2.42)), mse = c(29.44296, 32.2), cov = c(0.9440, 0.901), wid = c(20.14880, 19.0))



