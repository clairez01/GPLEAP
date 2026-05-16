# ------------------------------------------------------------------------------
# Program name: gpleap_sims.R
# Programmer..: Claire Zhu
# Date created: Feb 07 2025
# Description.: Conduct simulation study on the GPLEAP method using data 
#               generated based on the COMPACT trials.
# ------------------------------------------------------------------------------

# Load Libraries
library(nimble)
library(tidyverse)
library(posterior)
library(bayesplot)
library(haven)
library(ggplot2)

## Set directory paths
save.dir <- "/work/users/c/l/clairez1/Paper2/GPLEAP"
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
# ------------------------------------------------------------------------------
# Setup code for GPLEAP Model
# ------------------------------------------------------------------------------

# --------------------------------------------------
# Custom density for sampling f[1:n_total] | z
# --------------------------------------------------
dmix_mnorm <- nimbleFunction(
  run = function(x = double(1), z = double(1), Dsq_data = double(2),
                 alpha1 = double(0), rho1 = double(0),
                 alpha2 = double(0), rho2 = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    
    # Subset data based on z
    f1 <- x[z == 1]
    f2 <- x[z == 2]
    n1 <- sum(z[z == 1])
    n2 <- length(z) - n1
    Dsq_data1 <- Dsq_data[z == 1, z == 1]
    Dsq_data2 <- Dsq_data[z == 2, z == 2]
    
    # Compute parameters
    alpha1_sq     <- alpha1^2
    rho1_sq       <- rho1^2
    log_alpha1_sq <- log(alpha1_sq)
    alpha2_sq     <- alpha2^2
    rho2_sq       <- rho2^2
    log_alpha2_sq <- log(alpha2_sq)
    
    logProb <- 0
    
    # Compute C1
    if (n1 > 0) {
      C1 <- matrix(0, n1, n1)
      if (n1 > 1) {
        for (i in 1:(n1 - 1)) {
          C1[i, i] <- alpha1_sq + 1e-4
          for (j in (i + 1):n1) {
            C1[i, j] <- exp(log_alpha1_sq - Dsq_data1[i, j] / (2 * rho1_sq))
            C1[j, i] <- C1[i, j]
          }
        }
        C1[n1, n1] <- alpha1_sq + 1e-4
      } else if (n1 == 1) {
        C1[1, 1] <- alpha1_sq + 1e-4
      }
      mu1      <- rep(0, n1)
      chol_C1  <- chol(C1)
      logProb1 <- dmnorm_chol(x = f1, mean = mu1, cholesky = chol_C1, prec_param = FALSE, log = TRUE)
      logProb  <- logProb + logProb1
    }
    
    # Compute C2
    if (n2 > 0) {
      C2 <- matrix(0, n2, n2)
      if (n2 > 1) {
        for (i in 1:(n2 - 1)) {
          C2[i, i] <- alpha2_sq + 1e-4
          for (j in (i + 1):n2) {
            C2[i, j] <- exp(log_alpha2_sq - Dsq_data2[i, j] / (2 * rho2_sq))
            C2[j, i] <- C2[i, j]
          }
        }
        C2[n2, n2] <- alpha2_sq + 1e-4
      } else if (n2 == 1) {
        C2[1, 1] <- alpha2_sq + 1e-4
      }
      mu2      <- rep(0, n2)
      chol_C2  <- chol(C2)
      logProb2 <- dmnorm_chol(x = f2, mean = mu2, cholesky = chol_C2, prec_param = FALSE, log = TRUE)
      logProb  <- logProb + logProb2
    }
    
    if (log) return(logProb)
    else return(exp(logProb))
  })

# --------------------------------------------------
# Model code
# --------------------------------------------------
gpleap_code <- nimbleCode({
  # Priors for hyperparameters
  gamma         ~ dbeta(p_shape1, p_shape2) # exchangeability parameter
  log(rho1)     ~ dnorm(0, sd = 1)   
  log(alpha1)   ~ dnorm(0, sd = 1)     
  log(sigmasq1) ~ dnorm(0, sd = 1)
  log(rho2)     ~ dnorm(0, sd = 1)   
  log(alpha2)   ~ dnorm(0, sd = 1)     
  log(sigmasq2) ~ dnorm(0, sd = 1)
  omega[1:n]    ~ ddirch(ones[1:n]) # for bayesian bootstrap
  
  # Latent exchangeability indicators (for historical data)
  zprob[1:2] <- c(gamma, 1 - gamma)
  for(i in 1:n0) {
    z0[i] ~ dcat(prob = zprob[1:2])
  }
  # Current data are always exchangeable
  z[1:n_total] <- c(ones[1:n], z0[1:n0])
  
  # Distribution of f
  f[1:n_total] ~ dmix_mnorm(z = z[1:n_total], Dsq_data = Dsq_data[1:n_total, 1:n_total], 
                            alpha1 = alpha1, rho1 = rho1, alpha2 = alpha2, rho2 = rho2)
  
  # Likelihood
  for (i in 1:n_total) {
    z_bin[i]   <- 2 - z[i] # Turn z into a vector of 1s and 0s
    sigmasq[i] <- z_bin[i] * sigmasq1 + (1 - z_bin[i]) * sigmasq2
    y[i] ~ dnorm(f[i], var = sigmasq[i])
  }
  
  # Placeholders 
  # (needed for custom sampler but not directly used in the model code)
  placeholder1[1:N, 1:N]             <- Dsq_full[1:N, 1:N]
  placeholder2[1:n_total, 1:n_total] <- Dsq_data[1:n_total, 1:n_total]
  y_curr[1:n]                        <- curr_y[1:n]
})

# --------------------------------------------------
# Custom sampler for f | Y
# --------------------------------------------------
f_sampler_gpleap <- nimbleFunction(
  name = 'f_sampler_gpleap'  ## give name for sampler
  , contains = sampler_BASE  ## doesn't change
  , setup = function(model, mvSaved, target, control) {  ## arguments always the same
    ## Place only things that are fixed here
  }
  ## This is the function that runs when the sampler is called
  , run = function() {
    y        <- model[['y']]
    Dsq_data  <- model[['Dsq_data']]
    alpha1   <- model[['alpha1']]
    rho1     <- model[['rho1']]
    sigmasq1 <- model[['sigmasq1']]
    
    n_total <- length(y)
    n       <- length(model$y_curr)
    z       <- model[['z']]
    n1      <- sum(z == 1)
    n2      <- n_total - n1
    f_new   <- numeric(n_total)
    
    # Compute C1
    y1            <- y[z == 1]
    Dsq_data1      <- Dsq_data[z == 1, z == 1]
    I_n1          <- diag(n1)
    alpha1_sq     <- alpha1^2
    rho1_sq       <- rho1^2
    log_alpha1_sq <- log(alpha1^2)
    
    C1 <- matrix(0, n1, n1) # Initialize C1 matrix
    
    for (i in 1:(n1 - 1)) {
      C1[i,i] <- alpha1_sq + 1e-4
      for (j in (i + 1):n1) {
        C1[i,j] <- exp( log_alpha1_sq - Dsq_data1[i,j] / (2 * rho1_sq))
        C1[j,i] <- C1[i,j]
      }
    }
    C1[n1,n1] <- alpha1_sq + 1e-4
    
    y_prec1      <- C1 + sigmasq1 * I_n1
    y_Uprec1     <- chol(y_prec1)
    y_Uprec_inv1 <- backsolve(y_Uprec1, I_n1)
    y_cov1       <- y_Uprec_inv1 %*% t(y_Uprec_inv1)
    mu1          <- (C1 %*% y_cov1 %*% y1)[,1]
    cov1         <- C1 - C1 %*% y_cov1 %*% C1
    chol_cov1    <- chol(cov1)
    
    ## sample f1 | Y
    f1 <- rmnorm_chol(1, mean = mu1, cholesky = chol_cov1, prec_param = FALSE)
    f_new[z == 1] <- f1
    
    if (n2 > 0) {
      y2            <- y[z == 2]
      Dsq_data2      <- Dsq_data[z == 2, z == 2]
      I_n2          <- diag(n2)
      alpha2        <- model[['alpha2']]
      rho2          <- model[['rho2']]
      sigmasq2      <- model[['sigmasq2']]
      alpha2_sq     <- alpha2^2
      rho2_sq       <- rho2^2
      log_alpha2_sq <- log(alpha2^2)
      
      C2 <- matrix(0, n2, n2) # Initialize C1 matrix
      
      if (n2 > 1) {
        for ( i in 1:(n2 - 1) ) {
          C2[i,i] <- alpha2_sq + 1e-4
          for ( j in (i+1):n2 ) {
            C2[i,j] <- exp( log_alpha2_sq - Dsq_data2[i,j] / (2 * rho2_sq))
            C2[j,i] <- C2[i,j]
          }
        }
        C2[n2,n2] <- alpha2_sq + 1e-4
      }
      else {
        C2[n2,n2] <- alpha2_sq + 1e-4
      }
      
      y_prec2      <- C2 + sigmasq2 * I_n2
      y_Uprec2     <- chol(y_prec2)
      y_Uprec_inv2 <- backsolve(y_Uprec2, I_n2)
      y_cov2       <- y_Uprec_inv2 %*% t(y_Uprec_inv2)
      mu2          <- (C2 %*% y_cov2 %*% y2)[,1]
      cov2         <- C2 - C2 %*% y_cov2 %*% C2
      chol_cov2    <- chol(cov2)
      
      ## sample f2 | Y
      f2 <- rmnorm_chol(1, mean = mu2, cholesky = chol_cov2, prec_param = FALSE)
      f_new[z == 2] <- f2
    }
    
    ## Store sampled value and recalculate (note: <<- is necessary to store samples)
    model[['f']] <<- f_new
    model$calculate(target)  ## always needed at end
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)  ## always needed at end
  }
  ## Leave this as is
  , methods = list(
    reset = function() {}
  )
)

# --------------------------------------------------
# Compute treatment effect with compiled function
# --------------------------------------------------
trtEffect_GPLEAP_Batch <- nimbleFunction(
  setup = function(model, param_names_list) {
    # Move param_names into setup to lock them for compilation
    p_names <- param_names_list
    np <- length(p_names)
    
    # Pre-calculate dimensions
    N_total <- dim(model$Dsq_full)[1]
    n_curr  <- length(model$y_curr)
    n_tot   <- length(model$y)
    n_pred  <- N_total - n_tot
    Q_val   <- as.integer((n_pred / n_curr))
  },
  
  run = function(samples = double(2)) {
    n_samples <- dim(samples)[1]
    Q         <- Q_val
    n         <- length(model$y_curr)
    n_total   <- length(model$y)
    N         <- dim(model$Dsq_full)[1]
    
    # Storage matrix
    TE_matrix <- matrix(0, nrow = n_samples, ncol = Q)
    for (j in 1:n_samples) {
      values(model, p_names) <<- samples[j, 1:np]
      model$calculate() 
      
      # Extract current model values
      alpha1   <- model[['alpha1']]
      rho1     <- model[['rho1']]
      z        <- model[['z']]
      f        <- model[['f']]
      Dsq_full <- model[['Dsq_full']]
      omega    <- model[['omega']]
      
      # GP Logic
      n1        <- sum(z == 1)
      Z         <- c(z, rep(1, (N - n_total)))
      Dsq_full1 <- Dsq_full[Z == 1, Z == 1]
      f1        <- f[z == 1]
      K_size <- dim(Dsq_full1)[1]
      K      <- matrix(0, nrow = K_size, ncol = K_size)
      
      for(i in 1:K_size) {
        K[i,i] <- alpha1^2 + 0.001      
        for(k in 1:(i-1)) {
          K[i,k] <- alpha1^2 * exp(- Dsq_full1[i,k] / (2*rho1^2) )
          K[k,i] <- K[i,k]
        }
      }
      
      # Conditional Sampling
      U        <- chol(K)
      m        <- K_size - n1
      w_solve  <- forwardsolve(t(U[1:n1, 1:n1]), f1[1:n1])
      mu_fstar <- t(U[1:n1, (n1+1):K_size]) %*% w_solve[1:n1]
      zstar    <- rnorm(m)
      fstar    <- mu_fstar + U[(n1+1):K_size, (n1+1):K_size] %*% zstar
      
      # 5. Bayesian-bootstrap marginal means
      mu <- numeric(Q)
      for (q in 1:Q) {
        mu[q] <- inprod(omega[1:n], fstar[(n*q - n + 1):(n*q), 1])
      }
      
      # Compute treatment effects
      sum_treated <- 0
      for (q in 1:(Q-1)) {
        TE_matrix[j, q] <- mu[q] - mu[Q] 
        sum_treated     <- sum_treated + mu[q]
      }
      
      avg_treated <- sum_treated / (Q - 1)
      TE_matrix[j, Q] <- avg_treated - mu[Q] 
    }
    
    returnType(double(2))
    return(TE_matrix)
  }
)

# End setup code

################################################################################

current_i <- subset(current, iteration == 1)
historical_i <- subset(historical, q == q.id & prob.exch == prob.exch.id & setting == data.type & n0 == n0.id)

# Assign data to y, X, y0, and X0
n  <- 80
y  <- switch(data.type.id, current_i$y1, current_i$y2, current_i$y3)
X  <- cbind(current_i$dose)
y0 <- historical_i$y0
X0 <- cbind(historical_i$dose)

sq_euclidean <- function(A, B) {
  A_sq <- rowSums(A^2)
  B_sq <- rowSums(B^2)
  outer(A_sq, B_sq, "+") - 2 * A %*% t(B)
}

# ------------------------------
# Assign initial data values to compile the models
# continuous dose
Q       <- 3
n       <- length(y)
n0      <- length(y0)
n_total <- n + n0
N       <- n_total + Q*n

X_all   <- rbind(X, X0)
X_high  <- as.matrix(current_i$weight*60/1000)
X_low   <- as.matrix(current_i$weight*40/1000)
X_ctrl  <- as.matrix(rep(0, n))

X_full   <- rbind(X_all, X_high, X_low, X_ctrl)
Dsq_full <- sq_euclidean(X_full, X_full)
Dsq_data <- sq_euclidean(X_all, X_all)

# --------------------------
# Compile GPLEAP Model
# --------------------------

gpleap_data   <- list(y = c(y, y0), curr_y = y, Dsq_full = Dsq_full, Dsq_data = Dsq_data
                      , ones = rep(1, n), p_shape1 = 1, p_shape2 = 1)
gpleap_consts <- list(n = n, n0 = n0, n_total = n + n0, N = n_total + Q*n)
gpleap_inits  <- list(sigmasq1 = 1, sigmasq2 = 1, f = rnorm(n_total), omega = rep(1/n, n)
                      , z0 = rbinom(n0, 1, 0.5) + 1)

# Create / configure model
gpleap_model  <- nimbleModel(gpleap_code, data = gpleap_data, inits = gpleap_inits, constants = gpleap_consts)
gpleap_cmodel <- compileNimble(gpleap_model)
gpleap_conf   <- configureMCMC(gpleap_cmodel, monitors = c('gamma', 'sigmasq1', 'rho1', 'alpha1','f', 'omega', 'z'), print = FALSE)

## Change samplers
gpleap_conf$removeSampler('f')
gpleap_conf$addSampler('f', 'f_sampler_gpleap')

# Build model
gpleap_mcmc   <- buildMCMC(gpleap_conf)
gpleap_cmcmc  <- compileNimble(gpleap_mcmc, project = gpleap_model)

# Compile function for treatment effect
gpleap.test <- runMCMC(gpleap_cmcmc, niter = 1)
param_names_gpleap <- colnames(gpleap.test)

param_names_list <- colnames(gpleap.test)
batch_gpleap <- trtEffect_GPLEAP_Batch(gpleap_model, param_names_list)
C_batch_gpleap  <- compileNimble(batch_gpleap, project = gpleap_model)
te_names <- c("TE_HighDose", "TE_LowDose", "TE_Pooled")

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
  X_all  <- as.matrix(c(current_i$dose, historical_i$dose))
  y      <- switch(data.type.id, current_i$y1, current_i$y2, current_i$y3)
  X_high <- as.matrix(current_i$weight*60/1000)
  X_low  <- as.matrix(current_i$weight*40/1000)
  
  X_full   <- rbind(X_all, X_high, X_low, X_ctrl)
  Dsq_full <- sq_euclidean(X_full, X_full)
  Dsq_data <- sq_euclidean(X_all, X_all)
  
  # Replace data in nimble models
  gpleap_cmodel$y        <- c(y, y0)
  gpleap_cmodel$Dsq_data <- Dsq_data
  gpleap_cmodel$Dsq_full <- Dsq_full
  
  # WRAPPER START: Try to run MCMC and summarize results
  iter_res <- tryCatch({
    
    gpleap.smpl <- runMCMC(gpleap_cmcmc, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
    
    gpleap.summary <- data.frame(summarise_draws(gpleap.smpl, "mean", "median", "sd", "mad",
                                                 q2.5 = ~quantile(.x, probs = 0.025),  q97.5 = ~quantile(.x, probs = 0.975), 
                                                 "rhat", "mcse_mean",  "ess_bulk", "ess_tail" ))
    
    trt_effect_gpleap <- C_batch_gpleap$run(gpleap.smpl)
    te_draws <- as_draws_df(trt_effect_gpleap)
    variables(te_draws) <- te_names
    
    gpleap.trteffect  <- data.frame(summarise_draws(te_draws, "mean", "median", "sd", "mad",
                                                    q2.5 = ~quantile(.x, probs = 0.025),  q97.5 = ~quantile(.x, probs = 0.975), 
                                                    "rhat", "mcse_mean",  "ess_bulk", "ess_tail" ))
    
    # Combine results
    temp_res <- rbind(gpleap.summary[gpleap.summary$variable == "gamma", ], gpleap.trteffect)
    temp_res$iteration_id <- i  # Good for tracking which sims worked
    temp_res$status <- "success"
    
    temp_res # This is what gets assigned to iter_res
    
  }, error = function(e) {
    # What to do if NIMBLE crashes
    message(paste("Error in Iteration", i, ":", e$message))
    
    # Create a dummy row to record the failure
    fail_res <- data.frame(variable = "FAILED", mean = NA, iteration_id = i, status = "failed")
    return(fail_res)
  })
  # WRAPPER END
  
  # Append results if it wasn't a total failure, or keep the failure record
  iter_res$method <- "gpleap"
  results.all <- bind_rows(results.all, iter_res)
  
  cat("Iteration:", i, "| Status:", unique(iter_res$status), "\n")
}

end <- Sys.time()
print(end - start)

# Optional: Check how many failed at the end
# sum(results.all$variable == "FAILED")

# End sim code
# ------------------------------------------

## SAVE THE RESULTS

lst <- list(
  'simscen' = grid.id
  , 'id'      = id
  , 'simres'  = results.all
)

saveRDS(lst, filename)

