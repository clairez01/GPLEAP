## Directory to store files
save.dir <- '/work/users/c/l/clairez1/Paper2'

# --- Simulation Parameters ---
n0_list     <- c(15, 20, 30, 40, 50)  # Historical sample sizes
q_bias      <- c(-18, -9, 9, 18)      # Vertical shifts for non-exchangeable data
prob_exch   <- c(0, 0.2, 0.5, 0.8)    # Partial exchangeability levels
data_gen    <- c("linear", "quadratic", "emax") # Functional truths

# --- Create Unified Grid ---

# Grid 1: Fully Exchangeable (p = 1)
# In this case, q is effectively 0 because the trials match perfectly.
grid1 <- expand.grid(
  q         = 0, 
  prob.exch = 1, 
  data.type = data_gen,
  n0        = n0_list
)

# Grid 2: Partial/No Exchangeability (p < 1)
# Here we test all combinations of p, q, n0, and data.type
grid2 <- expand.grid(
  q         = q_bias, 
  prob.exch = prob_exch, 
  data.type = data_gen,
  n0        = n0_list
)

# Combine and Sort
grid <- rbind(grid1, grid2)
# Sorting by data type, then prob.exch, then bias, then sample size
grid <- grid[order(grid$data.type, -grid$prob.exch, grid$q, grid$n0), ]

# --- Expansion for Parallel Processing ---
# Create 20 rows of each scenario to divide up sims for cluster arrays
grid_expanded <- grid[rep(1:nrow(grid), each = 10), ]

# Unique Seed for each divided-up scenario
set.seed(4030)
grid_expanded$seed <- sample(seq_len(1000 * nrow(grid_expanded)), nrow(grid_expanded), replace = FALSE)

# Save the final grid
saveRDS(grid_expanded, file = file.path(save.dir, "grid_updated.rds"))
