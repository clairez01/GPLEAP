# ------------------------------------------------------------------------------
# Program name: comparison_sims.R
# Programmer..: Claire Zhu
# Date created: 10 Feb 2026
# Description.: Conduct simulation study on data generated based on the COMPACT 
#               trials using the PP and LEAP
# ------------------------------------------------------------------------------

# Load Libraries
library(tidyverse)
library(posterior)
library(bayesplot)
library(haven)
library(ggplot2)
library(cmdstanr)
library(hdbayes)

## Set directory paths
save.dir <- "/work/users/c/l/clairez1/Paper2/LEAP_PP"
data.dir <- "/work/users/c/l/clairez1/Paper2"

# Load in simulation data
setwd(data.dir)
current <- readRDS("simdata_current.rds")
historical <- readRDS("simdata_historical.rds")

# Load in grid of simulation parameters
grid <- readRDS("/work/users/c/l/clairez1/Paper2/grid_updated.rds")
#grid <- subset(grid, n0 == 40)
grid <- subset(grid, n0 != 40)

## get task ID
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if ( is.na(id ) )
  id <- 1

## get this job's sim parameters
grid.id      <- grid[id, ]
seed.id      <- as.integer(grid.id$seed)
n0.id        <- as.numeric(grid.id$n0)
q.id         <- as.numeric(grid.id$q)
prob.exch.id <- as.numeric(grid.id$prob.exch)
data.type.id <- as.numeric(grid.id$data.type)
data.type    <- grid.id$data.type

## Obtain file name based on id
filename <- file.path(save.dir, paste0('id_', id, '_n0_', n0.id, '_q_', q.id, 
                                       '_probexch_', prob.exch.id, '_datatype_', data.type.id, '.rds'))

set.seed(seed.id)

################################################################################
current_i <- subset(current, iteration == 1)
historical_i <- subset(historical, q == q.id & prob.exch == prob.exch.id & setting == data.type & n0 == n0.id)

# Assign data to y, X, y0, and X0
n  <- nrow(current_i)
n0 <- nrow(historical_i)
y  <- switch(data.type.id, current_i$y1, current_i$y2, current_i$y3)
X  <- cbind(current_i$dose)
y0 <- historical_i$y0
X0 <- historical_i$dose #rep(1, n0)
df0 <- data.frame(outcome = as.numeric(y0), V2 = as.numeric(X0))

calc_leap_linear_ate <- function(leap_smpl, q_target, w_current, scale_factor = 1000) {
  
  # 1. Extract the posterior draws for the dose coefficient (slope)
  beta1_draws <- leap_smpl$V2
  n_draws <- length(beta1_draws)
  n_current <- length(w_current)
  
  # 2. Calculate the scaled continuous target doses for all current subjects
  v2_target <- (q_target * w_current) / scale_factor
  
  # Initialize vector to store the marginalized ATE for each MCMC draw
  ate_draws <- numeric(n_draws)
  
  # 3. Apply the Bayesian Bootstrap marginalization
  for (m in 1:n_draws) {
    
    # Draw Bayesian Bootstrap weights (Dirichlet(1_n)) via normalized exponentials
    bb_weights <- rexp(n_current)
    bb_weights <- bb_weights / sum(bb_weights)
    
    # Marginalize over the empirical weight distribution
    ate_draws[m] <- beta1_draws[m] * sum(bb_weights * v2_target)
  }
  
  # 4. Reconstruct the draws format using the original MCMC metadata
  # This preserves the chain and iteration structures so diagnostics work properly
  draws_df_obj <- data.frame(
    .chain     = leap_smpl$.chain,
    .iteration = leap_smpl$.iteration,
    .draw      = leap_smpl$.draw,
    ATE        = ate_draws
  )
  draws_obj <- as_draws_df(draws_df_obj)
  
  # 5. Format the output as a dataframe row matching leap.summary
  summary_row <- data.frame(
    variable  = paste0("ATE_q", q_target),
    mean      = mean(ate_draws),
    median    = median(ate_draws),
    sd        = sd(ate_draws),
    mad       = mad(ate_draws),
    `2.5%`    = unname(quantile(ate_draws, probs = 0.025)),
    `97.5%`   = unname(quantile(ate_draws, probs = 0.975)),
    rhat      = tryCatch(unname(rhat(draws_obj)), error = function(e) NA),
    mcse_mean = tryCatch(unname(mcse_mean(draws_obj)), error = function(e) NA),
    ess_bulk  = tryCatch(unname(ess_bulk(draws_obj)), error = function(e) NA),
    ess_tail  = tryCatch(unname(ess_tail(draws_obj)), error = function(e) NA),
    check.names = FALSE # Prevents R from changing '2.5%' to 'X2.5.'
  )
  
  return(summary_row)
}

# ------------------------------------------------------------------------------
# Begin simulation code

Nsims  <- 1000
niter  <- 5000
burnin <- 2000
thin   <- 1

# Create empty data frame to store simulation results
results.all <- data.frame()

start <- Sys.time()
for (i in 1:Nsims) {
  
  dfn <- i + Nsims * (ifelse(id %% 10 == 0, 10, id %% 10) - 1)
  current_i <- subset(current, iteration == dfn)
  
  # Assign new current data
  X  <- current_i$dose
  y  <- switch(data.type.id, current_i$y1, current_i$y2, current_i$y3)
  
  df <- data.frame(outcome = as.numeric(y), V2 = as.numeric(X)  )
  
  # --- LEAP MODEL ---
  leap.smpl <- glm.leap(formula = outcome ~ V2, family = gaussian(link = "identity"),
                        data.list = list(df, df0), K = 2, iter_warmup = burnin, iter_sampling = niter, chains = 1, refresh = 0)
  
  leap_ate <- calc_leap_linear_ate(leap_smpl = leap.smpl, q_target = 40, w_current = current_i$weight)
  
  leap.summary <- leap.smpl %>% rename_variables(trteffect = V2) %>%  subset_draws(variable = c("trteffect", "probs[1]")) %>%
    summarize_draws("mean", "median", "sd", "mad", ~quantile(.x, probs = c(0.025, 0.975)), "rhat", "mcse_mean", "ess_bulk", "ess_tail") %>%
    as.data.frame()
  
  leap.summary <- rbind(leap.summary, leap_ate)
  
  # --- PP MODEL (a0 = 0.2) ---
  pp_0.2.smpl <- glm.pp(formula = outcome ~ V2, family = gaussian(link = "identity"), data.list = list(df, df0), a0.vals = 0.2,
                        iter_warmup = burnin, iter_sampling = niter, chains = 1, refresh = 0)
  
  pp_0.2_ate <- calc_leap_linear_ate(leap_smpl = pp_0.2.smpl, q_target = 40, w_current = current_i$weight)
  
  pp_0.2.summary <- pp_0.2.smpl %>% rename_variables(trteffect = V2) %>%  subset_draws(variable = "trteffect") %>%
    summarize_draws("mean", "median", "sd", "mad", ~quantile(.x, probs = c(0.025, 0.975)), "rhat", "mcse_mean", "ess_bulk", "ess_tail") %>%
    as.data.frame()
  
  pp_0.2.summary <- rbind(pp_0.2.summary, pp_0.2_ate)
  
  # --- PP MODEL (a0 = 0.5) ---
  pp_0.5.smpl <- glm.pp(formula = outcome ~ V2, family = gaussian(link = "identity"), data.list = list(df, df0), a0.vals = 0.5,
                        iter_warmup = burnin, iter_sampling = niter, chains = 1, refresh = 0)
  
  pp_0.5_ate <- calc_leap_linear_ate(leap_smpl = pp_0.5.smpl, q_target = 40, w_current = current_i$weight)
  
  pp_0.5.summary <- pp_0.5.smpl %>% rename_variables(trteffect = V2) %>%  subset_draws(variable = "trteffect") %>%
    summarize_draws("mean", "median", "sd", "mad", ~quantile(.x, probs = c(0.025, 0.975)), "rhat", "mcse_mean", "ess_bulk", "ess_tail") %>%
    as.data.frame()
  
  pp_0.5.summary <- rbind(pp_0.5.summary, pp_0.5_ate)
  
  # --- PP MODEL (a0 = 0.8) ---
  pp_0.8.smpl <- glm.pp(formula = outcome ~ V2, family = gaussian(link = "identity"), data.list = list(df, df0), a0.vals = 0.8,
                        iter_warmup = burnin, iter_sampling = niter, chains = 1, refresh = 0)
  
  pp_0.8_ate <- calc_leap_linear_ate(leap_smpl = pp_0.8.smpl, q_target = 40, w_current = current_i$weight)
  
  pp_0.8.summary <- pp_0.8.smpl %>% rename_variables(trteffect = V2) %>%  subset_draws(variable = "trteffect") %>%
    summarize_draws("mean", "median", "sd", "mad", ~quantile(.x, probs = c(0.025, 0.975)), "rhat", "mcse_mean", "ess_bulk", "ess_tail") %>%
    as.data.frame()
  
  pp_0.8.summary <- rbind(pp_0.8.summary, pp_0.8_ate)
  
  # --- COMBINE RESULTS ---
  leap.results <- leap.summary
  leap.results$method <- "leap"
  
  pp_0.2.results <- pp_0.2.summary
  pp_0.2.results$method <- "pp_0.2"
  
  pp_0.5.results <- pp_0.5.summary
  pp_0.5.results$method <- "pp_0.5"
  
  pp_0.8.results <- pp_0.8.summary
  pp_0.8.results$method <- "pp_0.8"
 
  # Bind all current iterations
  results <- rbind(leap.results, pp_0.2.results, pp_0.5.results, pp_0.8.results)
  
  # Append to master dataframe
  results.all <- rbind(results.all, results)
  
  cat("Iteration:", i, "\n")
}

end <- Sys.time()
end - start
# End sim code
# ------------------------------------------

## SAVE THE RESULTS

lst <- list(
  'simscen' = grid.id
  , 'id'      = id
  , 'simres'  = results.all
)

saveRDS(lst, filename)