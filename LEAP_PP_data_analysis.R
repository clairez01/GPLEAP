# ------------------------------------------------------------------------------
# Program name: LEAP_PP_data_analysis.R
# Programmer..: Claire Zhu
# Date created: 11 Feb 2026
# Description.: Conduct analysis on data from COMPACT Phase 2 and Phase 3 
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


# Load in CSL data
setwd("/nas/longleaf/home/clairez1/Paper2/CSL Datasets")
current <- read_sas("csl830_current.sas7bdat")
historical <- read_sas("csl830_external.sas7bdat")

current$DOSE_scaled    <- current$DOSE/1000
historical$DOSE_scaled <- historical$DOSE/1000

current$trtgrp2 <- ifelse(current$TRT01AN == 2, 2, ifelse(current$TRT01AN == 1, 1, 0))
current$trtgrp2 <- factor(current$trtgrp2)

X <- current$trtgrp2

# Analysis of actual data
y  <- current$CHG
X  <- current$DOSE_scaled #current$trtgroup
y0 <- historical$CHG
X0 <- historical$DOSE_scaled #historical$trtgroup

df  <- data.frame(outcome = as.numeric(y), V2 = as.numeric(X)  )
df0 <- data.frame(outcome = as.numeric(y0), V2 = as.numeric(X0))

set.seed(4002)
leap.smpl <- glm.leap(formula = outcome ~ V2, family = gaussian(link = "identity"),
                      data.list = list(df, df0), K = 2, iter_warmup = 2000, iter_sampling = 50000, chains = 4  )

leap.smpl %>% summarize_draws()

leap.summary <- leap.smpl %>%
  # Rename V2 while keeping probs[1] as is
  rename_variables(trteffect = V2) %>%
  # Select both variables of interest
  subset_draws(variable = c("trteffect", "probs[1]")) %>%
  # Apply all summary statistics
  summarize_draws(
    "mean", "median", "sd", "mad", 
    ~quantile(.x, probs = c(0.025, 0.975)), 
    "rhat", "mcse_mean", "ess_bulk", "ess_tail"
  ) %>%
  as.data.frame()

leap.summary


pp.smpl <- glm.pp(
  formula = outcome ~ V2,
  family = gaussian(link = "identity"),
  data.list = list(df, df0),
  a0.vals = 0.8,
  iter_warmup = 2000,
  iter_sampling = 50000,
  chains = 1
)

pp.summary <- pp.smpl %>%
  # Rename V2 to trteffect for consistency across your models
  rename_variables(trteffect = V2) %>%
  # Isolate the treatment effect
  subset_draws(variable = "trteffect") %>%
  # Calculate the requested statistics
  summarize_draws(
    "mean", "median", "sd", "mad",
    ~quantile(.x, probs = c(0.025, 0.975)),
    "rhat", "mcse_mean", "ess_bulk", "ess_tail"
  ) %>%
  as.data.frame()

pp.summary


df  <- data.frame(outcome = as.numeric(y), V2 = as.numeric(X)  )
df0 <- data.frame(outcome = as.numeric(y0), V2 = as.numeric(X0))

npp.smpl <- lm.npp(
  formula = outcome ~ V2,
  data.list = list(df, df0),
  sigmasq.shape = 2.1,
  sigmasq.scale = 1.1,
  a0.shape1 = 1,
  a0.shape2 = 1,
  iter_warmup = 2000,
  iter_sampling = 50000,
  chains = 1
)

# Step 1: Estimate the log normalizing constants
# Note: You only pass the historical data (df0) to this helper function
lognc_estimates <- glm.npp.lognc(
  formula = outcome ~ V2,
  family = gaussian(link = "identity"),
  histdata = df0,
  a0 = 0.8
)

# Step 2: Fit the model using the outputs from Step 1
npp.smpl <- glm.npp(
  formula = outcome ~ V2,
  family = gaussian(link = "identity"),
  data.list = list(df, df0),
  a0.lognc = lognc_estimates$a0,      # Passes the a0 grid
  lognc = lognc_estimates$lognc,      # Passes the estimated constants
  iter_warmup = burnin, 
  iter_sampling = niter, 
  chains = 1
)
