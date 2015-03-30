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

source("simul_prep_functions.R")



# simulate the data


nsite <- 100
nspec <- 2
nyear <- 5
nrep <- 28


# we can do it like this

sim_list <- gen_sim_list(nsite, nspec, nyear, nrep)

mats <- sim_matrices(sim_list)

# or even this

# make a list of the types of simulations we want to do
sims <- sim_all(nsite, nspec, nyear, nrep, add_NA = TRUE)


# start making larger scale simulations


nsite <- seq(10, 200, 10)
nspec = 2
nyear = 5
nrep = 28

big_sim <- vector(mode = "list", length = length(nsite))
for( i in 1:length(nsite)){
big_sim[[i]] <- sim_all(nsite = nsite[i], nspec, nyear, nrep)
} 
  










