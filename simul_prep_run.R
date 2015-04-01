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


nsite <- 118
nspec <- 5
nyear <- 10
nrep <- 28
gam <- c(.1, .5, .9)
phi <- c(.1, .5, .9)
p <- c(.1, .5, .9)
gam_sd <- c(.1, 1, 3)
phi_sd <- c(.1, 1, 3)
p_sd <- c(.1, 1, 3)

sim_df <- expand.grid(nsite, nspec, nyear, nrep, gam, phi, p, gam_sd, phi_sd, p_sd)
colnames(sim_df) <- c("nsite", "nspec", "nyear", "nrep", "gam",
                      "phi", "p", "gam_sd", "phi_sd", "p_sd")


# we can do it like this

sim_list <- gen_sim_list(nsite, nspec, nyear, nrep)

mats <- sim_matrices(sim_list)

# or even this

# make a list of the types of simulations we want to do
sims <- sim_all(nsite, nspec, nyear, nrep, add_NA = TRUE)


# start making larger scale simulations


nsite <- seq(10, 200, 10)
nspec = 10
nyear = 5
nrep = 28

#change the name of the portion of the list 'name'




    big_sim <- vector(mode = "list", length = nrow(sim_df))
    for( i in 1:nrow(sim_df)){
    big_sim[[i]] <- list(sim_all(nsite = sim_df$nsite[i], sim_df$nspec[i],
                            sim_df$nyear[i], sim_df$nrep[i],
                            actual_gam = sim_df$gam[i], sd_gam = sim_df$gam_sd[i],
                            actual_phi = sim_df$phi[i], sd_phi = sim_df$phi_sd[i],
                            actual_p = sim_df$p[i], sd_p = sim_df$p_sd[i]))
    names(big_sim[[i]]) <- paste("gam(", sim_df$gam[i],",", sim_df$gam_sd[i],
                                 ")_phi(", sim_df$phi[i], ",", sim_df$phi_sd[i],
                                 ")_p(", sim_df$p[i], ",", sim_df$p_sd[i], ")",
                                 sep = "")
                                 
                              
    }
    

sim_samp <- big_sim[sample(1:729, 2)]


  









