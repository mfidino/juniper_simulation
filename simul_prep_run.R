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

# load the functions to run everything
source("simul_prep_functions.R")



# simulate the data


# for our comparison stuff
sim_df <- make_sim_df(nspec = 3, nsite = 80, nyear = 6, nrep = 28,
                      gam = .4, phi = .7, p = .1, gam_sd = 1, beta = .5,
                      beta_sd = 1, beta2 = .5, beta_sd2 = .5,  inxs = 0, 
                      inxs_sd = 2, phi_sd = 1, p_sd = 1, add_NA = FALSE, percent_to_NA = 0.2,
                      cut_off = 20,
                      row_replicate = 2)



# test gens 2 simulations

#Note: takes ~0.5 seconds to run for each row of sim_df
#Note: known_inxs fills row wise, give number of inxs equal
#      to size of the entire inxs matrix (we write over the diagonal)
all_sim <- simulate_from_sim_df(sim_df, test = FALSE,
                                known_inxs = c(1,1,2,3,0,1,2,3,0))

# parameters to follow.
params <- c("gam",
            "phi",
            "psi_in", 
            "p_p",
            "area_gam",
            "area_phi",
            "dist2cent_gam",
            "dist2cent_phi",
            "p",
            "z")


  
start_time <- Sys.time()
batch_analyze(all_sim = list(all_sim[[1]]), params = params, n_chains = 7,adapt_steps = 3000,
                burn_in = 20000, sample_steps = 10000, thin_steps = 10,
              make_comparisons = TRUE, model = "interaction_model.R")
end_time <- Sys.time()
