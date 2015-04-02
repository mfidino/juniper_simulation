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
burn_in = 40000
sample_steps = 30000
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

for(i in 1:length(all_sim)){
  

  
  inits <- function(){ # Must be = instead of <-
    list( 
      z = all_sim[[i]]$mats$zinit,
      p_gam = rbeta(1,1,1),
      sigma_gam = runif(1, 0, 5),
      p_phi = rbeta(1,1,1),
      sigma_phi = runif(1, 0, 5),
      p_p = rbeta(1, 1, 1)
    )
  }
  

  
  mod_mcmc <- as.mcmc.list(run.jags( model="intercept_model.txt" , 
                                     monitor=params , 
                                     data=all_sim[[i]]$data_list ,  
                                     inits=inits , 
                                     n.chains=n_chains ,
                                     adapt=adapt_steps ,
                                     burnin=burn_in , 
                                     sample=ceiling(sample_steps / n_chains) ,
                                     thin=thin_steps ,
                                     summarise=FALSE ,
                                     plots=FALSE,
                                     method = "parallel"))

# writing out file stuff
basic_name <- base_file_name(all_sim[[i]])
  
write_diagnostics(mod_mcmc, iter = i, basic_name = basic_name)
write_summary(mod_mcmc, iter = i, basic_name = basic_name)


