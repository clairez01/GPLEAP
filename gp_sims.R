# ------------------------------------------------------------------------------
# Program name: gp_sims.R
# Programmer..: Claire Zhu
# Date created: 10 Feb 2026
# Description.: Conduct simulation study on data generated based on the COMPACT 
#               trials using the GP
# ------------------------------------------------------------------------------

# Load Libraries
library(nimble)
library(tidyverse)
library(posterior)
library(bayesplot)
library(haven)
library(ggplot2)

## Set directory paths
save.dir <- "/work/users/c/l/clairez1/Paper2/GP"
data.dir <- "/work/users/c/l/clairez1/Paper2"

# Load in simulation data
setwd(data.dir)
current <- readRDS("simdata_current.rds")

# Load in grid of simulation parameters
grid <- readRDS("/work/users/c/l/clairez1/Paper2/grid_updated.rds")
grid <- subset(grid, prob.exch == 1.0 & n0 != 40 )

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
# Setup code for GP model
# ------------------------------------------------------------------------------
gp_sample_conditional <- nimbleFunction(
  run = function( f = double(1), U = double(2) ) {
    n <- length(f)
    N <- dim(U)[1]
    m <- N - n
    
    # Compute conditional distribution for f_star
    w_solve  <- forwardsolve(t(U[1:n, 1:n]), f[1:n])   
    mu_fstar <- t( U[1:n, (n+1):N] ) %*% w_solve[1:n]
    zstar <- rnorm(m)
    fstar <- mu_fstar + U[(n+1):N, (n+1):N] %*% zstar
    
    ## Return sample from predictive density
    returnType(double(1))
    return(fstar[, 1])
  }
)

gp_code <- nimbleCode({
  # Priors
  log(rho)     ~ dnorm(0, sd = 1)   
  log(alpha)   ~ dnorm(0, sd = 1)     
  log(sigmasq) ~ dnorm(0, sd = 1)
  omega[1:n]   ~ ddirch(ones[1:n])

  ## Build joint squared-exponential kernel on stacked X_full 
  # X_full = (X_obs, X_dose1, X_dose2, X_control)
  # Dimension N = n + Qn where here Q = 2, so Qn = 2n, N = 3n.
  for(i in 1:N) {
    K[i,i] <- alpha^2 + 0.001      
    for(j in 1:(i-1)) {
      K[i,j] <- alpha^2 * exp(- Dsq_full[i,j] / (2*rho^2) )
      K[j,i] <- K[i,j]
    }
  }
  
  # Take the Cholesky
  U[1:N,1:N] <- chol(K[1:N,1:N])
  
  # GP prior for f 
  f[1:n] ~ dmnorm(zeros[1:n], cholesky = U[1:n, 1:n], prec_param = 0)
  
  # Likelihood 
  Sigma[1:n,1:n] <- sigmasq * eye[1:n, 1:n]
  y[1:n] ~ dmnorm(f[1:n], cov = Sigma[1:n,1:n])
  
  # Sample f* | f
  fstar[1:Qn] <- gp_sample_conditional(f[1:n], U[1:N, 1:N])
  
  # Bayesian-bootstrap marginal means
  for ( q in 1:Q ) {
    mu[q] <- inprod( omega[1:n], fstar[ (n*q - n + 1):(n*q) ] )
  }
  
  # Compute treatment effects
  for (q in 1:(Q-1)) {
    TE[q] <- mu[q] - mu[Q]
  }
  
  TE_pooled <- sum(TE[1:(Q-1)]) / (Q - 1)
})

################################################################################

current_i <- subset(current, iteration == 1)

# Assign data to y, X, y0, and X0
n  <- 80
y  <- switch(data.type.id, current_i$y1, current_i$y2, current_i$y3)
X  <- cbind(current_i$dose)

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
N       <- n + Q*n

X_high  <- as.matrix(current_i$weight*60/1000)
X_low   <- as.matrix(current_i$weight*40/1000)
X_ctrl  <- as.matrix(rep(0, n))

X_full   <- rbind(X, X_high, X_low, X_ctrl)
Dsq_full <- sq_euclidean(X_full, X_full)

# --------------------------
# Compile GP Model
# --------------------------

gp_data   <- list(y = y, Dsq_full = Dsq_full
                  , zeros = rep(0, n), ones = rep(1, n), eye = diag(1, n))
gp_consts <- list(n = n, Qn = Q*n, N = n + Q*n, Q = Q)
gp_inits  <- list(sigmasq = 1, f = rnorm(n), fstar = rep(mean(y), n*Q), omega = rep(1/n, n))

# Create / configure model
gp_model  <- nimbleModel(gp_code, data = gp_data, inits = gp_inits, constants = gp_consts)
gp_cmodel <- compileNimble(gp_model)
gp_conf   <- configureMCMC(gp_cmodel, monitors = c('sigmasq', 'rho', 'alpha', 'mu', 'TE', 'TE_pooled'), print = FALSE)

# Build model
gp_mcmc   <- buildMCMC(gp_conf)
gp_cmcmc  <- compileNimble(gp_mcmc, project = gp_model)

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
  X  <- cbind(current_i$dose)
  y  <- switch(data.type.id, current_i$y1, current_i$y2, current_i$y3)
  
  X_full   <- rbind(X, X_high, X_low, X_ctrl)
  Dsq_full <- sq_euclidean(X_full, X_full)
  
  # Replace data in nimble models
  gp_cmodel$y        <- y
  gp_cmodel$Dsq_full <- Dsq_full
  
  gp.smpl  <- runMCMC(gp_cmcmc, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
  gp.summary  <- data.frame(summarize_draws(gp.smpl, "mean", "median", "sd", "mad",
                ~quantile(.x, probs = c(0.025, 0.975)), "rhat", "mcse_mean",  "ess_bulk", "ess_tail")  )
  
  gp.results <- gp.summary[gp.summary$variable %in% c("TE[1]", "TE[2]", "TE_pooled"), ]
  gp.results$method <- "gp"
  
  results.all <- rbind(results.all, gp.results)
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