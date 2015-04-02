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


sim_df <- make_sim_df(nspec = 5, nsite = 118, nyear = 10, nrep = 28,
                      gam = c(.1, .5, .9), phi = c(.1, .5, .9),
                      p = c(.1, .5, .9), gam_sd = c(.1, 1, 3),
                      phi_sd = c(.1, 1, 3), p_sd = c(.1, 1, 3),
                      add_NA = TRUE, percent_to_NA = 0.2)



# test gens 2 simulations

all_sim <- simulate_from_sd_df(sim_df, test = TRUE)
    



# naming convention
  
names(big_sim[[i]]) <- paste("gam(", sim_df$gam[i],",", sim_df$gam_sd[i],
                             ")_phi(", sim_df$phi[i], ",", sim_df$phi_sd[i],
                             ")_p(", sim_df$p[i], ",", sim_df$p_sd[i], ")",
                             sep = "")    









