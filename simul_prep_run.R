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


sim_df <- make_sim_df(nspec = 5, nsite = 118, nyear = 5, nrep = 28,
                      gam = c(.1, .5, .9), phi = c(.1, .5, .9),
                      p = c(.1, .5, .9), gam_sd = c(.1, 1, 3),
                      phi_sd = c(.1, 1, 3), p_sd = c(.1, 1, 3),
                      add_NA = FALSE, percent_to_NA = 0.2, row_replicate = 1)
# for our comparison stuff
sim_df <- make_sim_df(nspec = 5, nsite = 118, nyear = 5, nrep = 28,
                      gam = .4, phi = .7, p = .1, gam_sd = 1,
                      phi_sd = 1, p_sd = 1, add_NA = FALSE, percent_to_NA = 0.2,
                      row_replicate = 10)



# test gens 2 simulations

all_sim <- simulate_from_sim_df(sim_df, test = TRUE)


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


  
start_time <- Sys.time()
batch_analyze(all_sim = list(all_sim[[1]]), params = params, n_chains = 4,adapt_steps = 3000,
                burn_in = 15000, sample_steps = 25000, thin_steps = 1,
              make_comparisons = TRUE)
end_time <- Sys.time()
