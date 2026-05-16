#-------------------------------------------------------------------------------
# Generate data sets for simulations
#-------------------------------------------------------------------------------

# Load in libraries
library(haven)
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

# Plot quadratic function from dose 0 to 9
curve(beta2[1] + beta2[2]*x + beta2[3]*x^2, from = 0, to = 9, 
      col = "red", lwd = 2, 
      main = "Dose Response", xlab = "Dose", ylab = "Response")

# Plot Emax function from dose 0 to 9
curve(beta3[1] + (beta3[2]*x)/(beta3[3] + x), from = 0, to = 9, 
      col = "red", lwd = 2, 
      main = "Dose Response", xlab = "Dose", ylab = "Response")


#-------------------------------------------------------------------------------
# Generate three sets of current data
# 1) Linear dose relationship
# 2) Quadratic dose relationship
# 3) Emax dose relationship

# --- Define "true" mean functions ---
mean_fun_y1 <- function(dose, beta) {beta[1] + beta[2] * dose}
mean_fun_y2 <- function(dose, beta) {beta[1] + beta[2] * dose + beta[3] * dose^2}
mean_fun_y3 <- function(dose, beta) {beta[1] + (beta[2] * dose) / (beta[3] + dose)}

# Set parameters
n_reps <- 10000
n_curr <- 80 

# Pre-allocate the data frame with the final number of rows
# Total rows = iterations * samples per iteration
total_rows <- n_reps * n_curr
sim_data <- data.frame(
  iteration = integer(total_rows),
  dose      = numeric(total_rows),
  trtgroup  = character(total_rows),
  weight    = numeric(total_rows),
  y1        = numeric(total_rows),
  y2        = numeric(total_rows),
  y3        = numeric(total_rows),
  stringsAsFactors = FALSE
)

set.seed(4001)
for (i in 1:n_reps) {
  
  # Calculate row indices for this specific iteration
  start_idx <- (i - 1) * n_curr + 1
  end_idx   <- i * n_curr
  
  # Sample with replacement
  curr_idx <- sample(nrow(current), n_curr, replace = TRUE)
  df_curr  <- current[curr_idx, ]
  d        <- df_curr$DOSE_scaled
  
  # Generate outcomes
  y1 <- rnorm(n_curr, mean_fun_y1(d, beta1), sd1)
  y2 <- rnorm(n_curr, mean_fun_y2(d, beta2), sd2)
  y3 <- rnorm(n_curr, mean_fun_y3(d, beta3), sd3)
  
  # Assign directly to the pre-allocated rows
  sim_data$iteration[start_idx:end_idx] <- i
  sim_data$dose[start_idx:end_idx]      <- d
  sim_data$trtgroup[start_idx:end_idx]  <- df_curr$trtgroup
  sim_data$weight[start_idx:end_idx]    <- df_curr$WEIGHT
  sim_data$y1[start_idx:end_idx]        <- y1
  sim_data$y2[start_idx:end_idx]        <- y2
  sim_data$y3[start_idx:end_idx]        <- y3
}

# Check structure
str(sim_data)
head(sim_data)

# Save datasets
setwd("/work/users/c/l/clairez1/Paper2")
#saveRDS(sim_data, file = "simdata_current.rds")
saveRDS(sim_data, file = "simdata_current_misspec.rds")

#-------------------------------------------------------------------------------
# Generate historical data
# --- Parameters ---
n0.id_list <- c(15, 20, 30, 40, 50) # Updated to a vector of sample sizes
prob.id <- c(0, 0.2, 0.5, 0.8, 1)  
q.id <- c(-18, -9, 9, 18)   
niter_search <- 10000         

set.seed(4022)
hist_data <- data.frame()

# Track total run time
start_time <- Sys.time()
cat("Starting historical data generation at:", format(start_time), "\n\n")

for (setting in 1:3) {
  setting_name <- switch(setting, "linear", "quadratic", "emax")
  
  # Assign true DGP parameters based on the current setting
  beta_truth <- switch(setting, beta1, beta2, beta3)
  sd_truth   <- switch(setting, sd1, sd2, sd3)
  mean_fun   <- switch(setting, mean_fun_y1, mean_fun_y2, mean_fun_y3)
  
  for (p in prob.id) {
    q_values <- if (p == 1) 0 else q.id
    
    for (q in q_values) {
      beta_non_exch <- beta_truth
      beta_non_exch[1] <- beta_non_exch[1] + q
      
      # Iterate over each historical sample size
      for (n0 in n0.id_list) {
        
        min_param_dist <- Inf
        best_df <- NULL
        
        for (i in 1:niter_search) {
          n_exch <- round(n0 * p)
          n_unexch <- n0 - n_exch
          
          # ---------------------------------------------------------
          # 1. Generate Exchangeable Data
          # ---------------------------------------------------------
          if (n_exch > 0) {
            idx_exch <- sample(nrow(historical), n_exch, replace=TRUE)
            d_exch   <- historical$DOSE_scaled[idx_exch]
            w_exch   <- historical$WEIGHT[idx_exch] 
            mu_exch  <- mean_fun(d_exch, beta_truth)
            y_exch   <- rnorm(n_exch, mu_exch, sd_truth)
          } else {
            d_exch <- w_exch <- y_exch <- numeric(0)
          }
          
          # ---------------------------------------------------------
          # 2. Generate Non-Exchangeable Data
          # ---------------------------------------------------------
          if (n_unexch > 0) {
            idx_unexch <- sample(nrow(historical), n_unexch, replace=TRUE)
            d_unexch   <- historical$DOSE_scaled[idx_unexch]
            w_unexch   <- historical$WEIGHT[idx_unexch] 
            mu_unexch  <- mean_fun(d_unexch, beta_non_exch)
            y_unexch   <- rnorm(n_unexch, mu_unexch, sd_truth)
          } else {
            d_unexch <- w_unexch <- y_unexch <- numeric(0)
          }
          
          # Combine temporarily for model fitting
          d_temp <- c(d_exch, d_unexch)
          y_temp <- c(y_exch, y_unexch)
          
          # ---------------------------------------------------------
          # 3. Dynamic Parameter Recovery Search
          # ---------------------------------------------------------
          fit <- NULL
          total_dist <- Inf
          
          if (setting == 1) { 
            # --- LINEAR SETTING ---
            fit <- tryCatch(lm(y_temp ~ d_temp), error = function(e) NULL)
            
            if (!is.null(fit)) {
              est_p1 <- coef(fit)["(Intercept)"]
              est_p2 <- coef(fit)["d_temp"]
              est_sigma <- summary(fit)$sigma
              
              exp_p1 <- (p * beta_truth[1]) + ((1 - p) * beta_non_exch[1])
              exp_p2 <- (p * beta_truth[2]) + ((1 - p) * beta_non_exch[2])
              
              dist_p1 <- (est_p1 - exp_p1)^2 / (exp_p1^2 + 1e-8)
              dist_p2 <- (est_p2 - exp_p2)^2 / (exp_p2^2 + 1e-8)
              dist_sigma <- (est_sigma - sd_truth)^2 / (sd_truth^2 + 1e-8)
              
              total_dist <- dist_p1 + dist_p2 + dist_sigma
            }
            
          } else if (setting == 2) { 
            # --- QUADRATIC SETTING ---
            fit <- tryCatch(lm(y_temp ~ d_temp + I(d_temp^2)), error = function(e) NULL)
            
            if (!is.null(fit)) {
              est_p1 <- coef(fit)["(Intercept)"]
              est_p2 <- coef(fit)["d_temp"]
              est_p3 <- coef(fit)["I(d_temp^2)"]
              est_sigma <- summary(fit)$sigma
              
              exp_p1 <- (p * beta_truth[1]) + ((1 - p) * beta_non_exch[1])
              exp_p2 <- (p * beta_truth[2]) + ((1 - p) * beta_non_exch[2])
              exp_p3 <- (p * beta_truth[3]) + ((1 - p) * beta_non_exch[3])
              
              dist_p1 <- (est_p1 - exp_p1)^2 / (exp_p1^2 + 1e-8)
              dist_p2 <- (est_p2 - exp_p2)^2 / (exp_p2^2 + 1e-8)
              dist_p3 <- (est_p3 - exp_p3)^2 / (exp_p3^2 + 1e-8)
              dist_sigma <- (est_sigma - sd_truth)^2 / (sd_truth^2 + 1e-8)
              
              total_dist <- dist_p1 + dist_p2 + dist_p3 + dist_sigma
            }
            
          } else if (setting == 3) { 
            # --- EMAX SETTING ---
            exp_p1 <- (p * beta_truth[1]) + ((1 - p) * beta_non_exch[1])
            exp_p2 <- (p * beta_truth[2]) + ((1 - p) * beta_non_exch[2])
            exp_p3 <- (p * beta_truth[3]) + ((1 - p) * beta_non_exch[3])
            
            # Start values guide the algorithm to ensure convergence 
            fit <- tryCatch(
              nls(y_temp ~ E0 + (Emax * d_temp) / (ED50 + d_temp),
                  start = list(E0 = exp_p1, Emax = exp_p2, ED50 = exp_p3)),
              error = function(e) NULL
            )
            
            if (!is.null(fit)) {
              est_p1 <- coef(fit)["E0"]
              est_p2 <- coef(fit)["Emax"]
              est_p3 <- coef(fit)["ED50"]
              est_sigma <- summary(fit)$sigma
              
              dist_p1 <- (est_p1 - exp_p1)^2 / (exp_p1^2 + 1e-8)
              dist_p2 <- (est_p2 - exp_p2)^2 / (exp_p2^2 + 1e-8)
              dist_p3 <- (est_p3 - exp_p3)^2 / (exp_p3^2 + 1e-8)
              dist_sigma <- (est_sigma - sd_truth)^2 / (sd_truth^2 + 1e-8)
              
              total_dist <- dist_p1 + dist_p2 + dist_p3 + dist_sigma
            }
          }
          
          # ---------------------------------------------------------
          # 4. Save the "Best" Representative Dataset
          # ---------------------------------------------------------
          # Added !is.na(total_dist) to safely ignore unidentifiable models
          if (!is.null(fit) && !is.na(total_dist) && total_dist < min_param_dist) {
            min_param_dist <- total_dist
            best_df <- data.frame(
              dose    = d_temp,
              weight  = c(w_exch, w_unexch), 
              y0      = y_temp,
              is_exch = c(rep(1, n_exch), rep(0, n_unexch))
            )
          }
        }
        
        # Tag the best dataset with its scenario labels
        best_df$prob.exch <- p
        best_df$q         <- q
        best_df$setting   <- setting_name
        best_df$n0        <- n0 
        
        hist_data <- rbind(hist_data, best_df)
        
        # Print progress to the console
        cat(sprintf("[%s] Generated: %-9s | n0 = %2d | prob = %3.1f | q = %3d\n", 
                    format(Sys.time(), "%H:%M:%S"), setting_name, n0, p, q))
        flush.console() 
      }
    }
  }
}

cat("\nFinished! Total run time:", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes.\n")

# --- Final Cleanup ---
# Ensure setting is a factor for easier plotting/analysis later
hist_data$setting <- factor(hist_data$setting, levels = c("linear", "quadratic", "emax"))


str(hist_data)
table(hist_data$setting)
table(hist_data$prob.exch, hist_data$q)

# Save datasets
setwd("/work/users/c/l/clairez1/Paper2")
saveRDS(hist_data, "simdata_historical.rds") 

