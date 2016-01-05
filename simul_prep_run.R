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
                      gam = .4, phi = .7, p = 0.15, gam_sd = 1, beta = .5,
                      beta_sd = 1, beta2 = .5, beta_sd2 = .5,  inxs = 0, 
                      inxs_sd = 2, phi_sd = 1, p_sd = 0.5, add_NA = FALSE, percent_to_NA = 0.2,
                      cut_off = 20,
                      row_replicate = 400)



# test gens 2 simulations

#Note: takes ~0.5 seconds to run for each row of sim_df
#Note: known_inxs fills col wise, give number of inxs equal
#      to size of the entire inxs matrix (we write over the diagonal)
all_sim <- simulate_from_sim_df(sim_df, test = FALSE,
                                known_inxs = c(1,-0.5,1,-1,1,1,-1.5, -2,1))

# parameters to follow.
params <- c("gam",
            "phi",
            "psi_in", 
            "beta",
            "beta2",
            "p",
            "z",
            "sigma_int")

model <- "interaction_model_less_hyper_gamma.R"
load.module("glm")
start_time <- Sys.time()
batch_analyze(all_sim = all_sim, params = params, n_chains = 7,adapt_steps = 3000,
                burn_in = 10000, sample_steps = 10000, thin_steps = 10,
              make_comparisons = TRUE, model = model, m_name = "pcp")
end_time <- Sys.time()

pcp_time <- end_time - start_time
all_sim <- simulate_from_sim_df(sim_df, test = FALSE,
                                known_inxs = c(1,-1,-2,-0.5,1,-1,0, -0.7,1))

start_time <- Sys.time()
batch_analyze(all_sim = all_sim, params = params, n_chains = 7,adapt_steps = 3000,
              burn_in = 10000, sample_steps = 10000, thin_steps = 10,
              make_comparisons = TRUE, model = model,m_name = "hc")
end_time <- Sys.time()

hc_time <- end_time - start_time
