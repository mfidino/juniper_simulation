#########################
#
#
# simul_prep_run for simulation work
#
# 3/24/2015
#
# written by M. Fidino


#packages needed

library(plyr)
library(reshape2)
library(rjags)
library(runjags)
library(mcmcplots)

source("simul_prep_functions.R")



# simulate the data


sim_df <- make_sim_df(nspec = 2, nsite = 50, nyear = 5, nrep = 28,
                      gam = c(.1, .5, .9), phi = c(.1, .5, .9),
                      p = c(.1, .5, .9), gam_sd = c(.1, 1, 3),
                      phi_sd = c(.1, 1, 3), p_sd = c(.1, 1, 3),
                      add_NA = TRUE, percent_to_NA = 0.2)



# test gens 2 simulations

all_sim <- simulate_from_sd_df(sim_df, test = TRUE)

n_chains = 5
adapt_steps = 5000
burn_in = 1000
sample_steps = 1000
thin_steps = 2

params <- c("gam",
            "phi",
            "psi_in",
            "p_gam",
            "sigma_gam",
            "p_phi", 
            "sigma_phi", 
            "p_p",
            "sigma_p",
            "mu_gam",
            "mu_phi",
            "p")


  

batch_analyze(all_sim = all_sim, params = params, n_chains = 5,adapt_steps = 1000,
                burn_in = 1000, sample_steps = 1000, thin_steps = 2)
