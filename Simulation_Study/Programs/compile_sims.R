# Load in libraries
library(haven)
library(tidyverse)
library(mvtnorm)
library(stats)

# Load in CSL data
setwd("/nas/longleaf/home/clairez1/Paper2/CSL Datasets")
current <- read_sas("csl830_current.sas7bdat")
historical <- read_sas("csl830_external.sas7bdat")

current$DOSE_scaled    <- current$DOSE/1000
historical$DOSE_scaled <- historical$DOSE/1000

# Analysis of actual data
y  <- current$CHG
dose <- current$DOSE_scaled
p  <- ncol(dose)

# --- Extract MLEs from Current Data ---

fit_lin <- lm(y ~ dose)
fit_quad <- lm(y ~ dose + I(dose^2))
fit_emax <- nls(CHG ~ E0 + (Emax * DOSE_scaled) / (ED50 + DOSE_scaled), 
                data = current, 
                start = list(E0 = coef(fit_lin)[1], Emax = max(current$CHG), ED50 = 4))

# "true" parameters
beta1 <- coef(fit_lin)
beta2 <- c(0, 28, -4) #coef(fit_quad)
beta3 <- c(-10, 100, 0.5) #coef(fit_emax)
sd1   <- sigma(fit_lin) 
sd2   <- 18.7
sd3   <- 18.7

# Load in simulation data
setwd("/work/users/c/l/clairez1/Paper2")
current <- readRDS("simdata_current_misspec.rds")
historical <- readRDS("simdata_historical_misspec_v2.rds")

# Define "true" mean functions
mean_fun_y1 <- function(dose, beta) {beta[1] + beta[2] * dose}
mean_fun_y2 <- function(dose, beta) {beta[1] + beta[2] * dose + beta[3] * dose^2}
mean_fun_y3 <- function(dose, beta) {beta[1] + (beta[2] * dose) / (beta[3] + dose)}

exch.mean1 <- c(mean(mean_fun_y1(60*current$weight/1000, beta1)), mean(mean_fun_y1(40*current$weight/1000, beta1))) - mean_fun_y1(0, beta1)
exch.mean2 <- c(mean(mean_fun_y2(60*current$weight/1000, beta2)), mean(mean_fun_y2(40*current$weight/1000, beta2))) - mean_fun_y2(0, beta2)
exch.mean3 <- c(mean(mean_fun_y3(60*current$weight/1000, beta3)), mean(mean_fun_y3(40*current$weight/1000, beta3))) - mean_fun_y3(0, beta3)
exch.mean  <- cbind(exch.mean1, exch.mean2, exch.mean3)


# Set directory for simulation data
setwd("/work/users/c/l/clairez1/Paper2")
# Load in grid of simulation parameters
grid <- readRDS("/work/users/c/l/clairez1/Paper2/grid_updated.rds")
#grid <- subset(grid, data.type == "linear")
#grid <- subset(grid, (data.type == "quadratic" | data.type == "emax") & n0 == 40)
grid <- subset(grid, n0 != 40)

setwd("/work/users/c/l/clairez1/Paper2/GPLEAP_sensitivity")
# List all .rds files in the directory
file_list <- list.files(pattern = "\\.rds$")
file_list <- file_list[file_list != "grid.rds"]

numbers <- as.numeric(gsub("id_(\\d+)_.*", "\\1", file_list))

# Sort file list based on the extracted numbers
sorted_file_list <- file_list[order(numbers, file_list)]

# Print the sorted file list
print(sorted_file_list)

sorted_file_list <- grep("n0_50", sorted_file_list, value = TRUE, invert = TRUE)

# ------------------------------------------------------------------------------
summarize_metric <- function(df, true_val, label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  # Ensure column names match your simulation output (e.g., X2.5. vs q5)
  bias    <- mean(df$mean - true_val, na.rm = TRUE)
  mse     <- mean((df$mean - true_val)^2, na.rm = TRUE)
  covprob <- mean(df$X2.5. <= true_val & df$X97.5. >= true_val, na.rm = TRUE)
  width   <- mean(df$X97.5. - df$X2.5., na.rm = TRUE)
  
  return(data.frame(
    mean = mean(df$mean), sd = mean(df$sd), 
    LB = mean(df$X2.5.), UB = mean(df$X97.5.), 
    bias = bias, MSE = mse, cov_prob = covprob, width = width, 
    method = label
  ))
}

# Run this for GPLEAP
combine_datasets <- function(start_id, end_id, grid, exch.mean) {
  dfs <- list() 
  
  for (id in start_id:end_id) {
    grid.id <- grid[id,]
    file <- paste0('id_', id, '_n0_', grid.id$n0, '_q_', grid.id$q, 
                   '_probexch_', grid.id$prob.exch, 
                   '_datatype_', as.numeric(grid.id$data.type), '.rds')
    
    if (file.exists(file)) {
      temp <- readRDS(file)
      
      # Filter to keep only successful runs
      success_data <- temp$simres %>% 
        dplyr::filter(status == "success")
      
      # Check if the file had zero successful runs. If so, skip it entirely.
      if (nrow(success_data) == 0) {
        message(sprintf("Skipping ID %d: All iterations failed.", id))
        next # Jumps immediately to the next 'id' in the loop
      }
      
      # Only add to the list if there are successful rows
      dfs[[as.character(id)]] <- success_data
      
      last_dt <- as.numeric(grid.id$data.type) 
    }
  }
  
  if (length(dfs) == 0) return(NULL)
  
  # 1. Combine all replicates for this scenario
  results.all <- dplyr::bind_rows(dfs)
  
  if (nrow(results.all) == 0) return(NULL)
  
  # 3. Set Truths based on the scenario's data.type
  true_h <- exch.mean[1, last_dt]
  true_l <- exch.mean[2, last_dt]
  true_p <- mean(exch.mean[, last_dt])
  
  # 4. Summarize each group
  sum_h <- summarize_metric(results.all %>% filter(variable == "TE_HighDose"), true_h, "GPLEAP_High")
  sum_l <- summarize_metric(results.all %>% filter(variable == "TE_LowDose"),  true_l, "GPLEAP_Low")
  sum_p <- summarize_metric(results.all %>% filter(variable == "TE_Pooled"),   true_p, "GPLEAP_Pooled")
  
  # Gamma summary
  g_df  <- results.all %>% filter(variable == "gamma")
  if(nrow(g_df) > 0) {
    sum_g <- data.frame(
      mean = mean(g_df$mean), sd = mean(g_df$sd), 
      LB = mean(g_df$X2.5.), UB = mean(g_df$X97.5.), 
      bias = NA, MSE = NA, cov_prob = NA, width = NA, method = "GPLEAP_gamma"
    )
  } else {
    sum_g <- NULL
  }
  
  # Combine into final summary table
  sim.summary <- dplyr::bind_rows(sum_h, sum_l, sum_p, sum_g)
  return(sim.summary)
}

# Run this for GP
combine_datasets <- function(start_id, end_id, grid, exch.mean) {
  dfs <- list() 
  
  for (id in start_id:end_id) {
    grid.id <- grid[id,]
    file <- paste0('id_', id, '_n0_', grid.id$n0, '_q_', grid.id$q, 
                   '_probexch_', grid.id$prob.exch, 
                   '_datatype_', as.numeric(grid.id$data.type), '.rds')
    
    if (file.exists(file)) {
      temp <- readRDS(file)
      dfs[[as.character(id)]] <- temp$simres
      last_dt <- as.numeric(grid.id$data.type) 
    }
  }
  
  if (length(dfs) == 0) return(NULL)
  
  # 1. Combine all replicates for this scenario
  results.all <- dplyr::bind_rows(dfs)
  
  if (nrow(results.all) == 0) return(NULL)
  
  # 3. Set Truths based on the scenario's data.type
  true_h <- exch.mean[1, last_dt]
  true_l <- exch.mean[2, last_dt]
  true_p <- mean(exch.mean[, last_dt])
  
  # 4. Summarize each group
  sum_h <- summarize_metric(results.all %>% filter(variable == "TE[1]"), true_h, "GP_High")
  sum_l <- summarize_metric(results.all %>% filter(variable == "TE[2]"),  true_l, "GP_Low")
  sum_p <- summarize_metric(results.all %>% filter(variable == "TE_Pooled"),   true_p, "GP_Pooled")
  
  # Gamma summary
  g_df  <- results.all %>% filter(variable == "gamma")
  if(nrow(g_df) > 0) {
    sum_g <- data.frame(
      mean = mean(g_df$mean), sd = mean(g_df$sd), 
      LB = mean(g_df$X2.5.), UB = mean(g_df$X97.5.), 
      bias = NA, MSE = NA, cov_prob = NA, width = NA, method = "GPLEAP_gamma"
    )
  } else {
    sum_g <- NULL
  }
  
  # Combine into final summary table
  sim.summary <- dplyr::bind_rows(sum_h, sum_l, sum_p, sum_g)
  return(sim.summary)
}

# Run this for LEAP and PP
summarize_metric <- function(df, true_val, label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  # Ensure column names match your simulation output (e.g., X2.5. vs q5)
  bias    <- mean(df$mean - true_val, na.rm = TRUE)
  mse     <- mean((df$mean - true_val)^2, na.rm = TRUE)
  covprob <- mean(df$`2.5%` <= true_val & df$`97.5%` >= true_val, na.rm = TRUE)
  width   <- mean(df$`97.5%` - df$`2.5%`, na.rm = TRUE)
  
  return(data.frame(
    mean = mean(df$mean), sd = mean(df$sd), 
    LB = mean(df$`2.5%`), UB = mean(df$`97.5%`), 
    bias = bias, MSE = mse, cov_prob = covprob, width = width, 
    method = label
  ))
}

combine_datasets <- function(start_id, end_id, grid, exch.mean) {
  dfs <- list() 
  
  for (id in start_id:end_id) {
    grid.id <- grid[id,]
    file <- paste0('id_', id, '_n0_', grid.id$n0, '_q_', grid.id$q, 
                   '_probexch_', grid.id$prob.exch, 
                   '_datatype_', as.numeric(grid.id$data.type), '.rds')
    
    if (file.exists(file)) {
      temp <- readRDS(file)
      dfs[[as.character(id)]] <- temp$simres
      last_dt <- as.numeric(grid.id$data.type) 
    }
  }
  
  if (length(dfs) == 0) return(NULL)
  
  # 1. Combine all replicates for this scenario
  results.all <- dplyr::bind_rows(dfs)
  
  if (nrow(results.all) == 0) return(NULL)
  
  # 3. Set Truths based on the scenario's data.type
  true_l <- exch.mean[2, last_dt]
  
  # 4. Summarize each group
  # LEAP
  sum_leap <- summarize_metric(results.all %>% filter(variable == "ATE_q40" & method == "leap"), true_l, "LEAP")
  
  # Power Prior
  sum_pp02 <- summarize_metric(results.all %>% filter(variable == "ATE_q40" & method == "pp_0.2"), true_l, "pp_0.2")
  sum_pp05 <- summarize_metric(results.all %>% filter(variable == "ATE_q40" & method == "pp_0.5"), true_l, "pp_0.5")
  sum_pp08 <- summarize_metric(results.all %>% filter(variable == "ATE_q40" & method == "pp_0.8"), true_l, "pp_0.8")

  # Gamma summary
  g_df  <- results.all %>% filter(variable == "probs[1]")
  if(nrow(g_df) > 0) {
    sum_g <- data.frame(
      mean = mean(g_df$mean), sd = mean(g_df$sd), 
      LB = mean(g_df$`2.5%`), UB = mean(g_df$`97.5%`), 
      bias = NA, MSE = NA, cov_prob = NA, width = NA, method = "LEAP_prob"
    )
  } else {
    sum_g <- NULL
  }
  
  # Combine into final summary table
  sim.summary <- dplyr::bind_rows(
    sum_leap, sum_g,
    sum_pp02, sum_pp05, sum_pp08)
  return(sim.summary)
}

# ------------------------------------------------------------------------------
compiled.results <- data.frame()
total_expected_ids <- 2040

for (start_id in seq(1, total_expected_ids, by = 10)) {
  
  # Set the range for this specific scenario
  end_id <- min(start_id + 9, total_expected_ids)
  
  # Run the summary function for this block of 20
  # Note: ensure grid and exch.mean are already loaded in your environment
  sim.summary <- combine_datasets(start_id, end_id, grid, exch.mean)
  
  # Only proceed if we actually found files in this range
  if (!is.null(sim.summary)) {
    
    # Extract the scenario parameters from the first ID of the block
    grid.id <- grid[start_id, ]
    
    # 3. Add identifying columns (Scenario Metadata)
    sim.summary$q        <- grid.id$q
    sim.summary$pr       <- grid.id$prob.exch
    sim.summary$dt       <- as.numeric(grid.id$data.type)
    sim.summary$n0       <- as.numeric(grid.id$n0)
    sim.summary$id_range <- paste0(start_id, "-", end_id)
    
    # 4. Append to the master data frame
    compiled.results <- rbind(compiled.results, sim.summary)
    
    # Status update for your console/log
    cat(sprintf("Scenario %s (IDs %d-%d): Processed successfully.\n", 
                paste0("q=", grid.id$q, "_pr=", grid.id$prob.exch), start_id, end_id))
  } else {
    cat(sprintf("Scenario IDs %d-%d: No files found. Skipping.\n", start_id, end_id))
  }
}


head(compiled.results)


setwd("/nas/longleaf/home/clairez1/Paper2/Results")
saveRDS(compiled.results, file = "compiled_results_gpleap_sens.rds")

