#########################
#
#
# simul_prep_run for simulation work
#
# 3/24/2015
#
# written by M. Fidino



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

sims <- sim_all(nsite, nspec, nyear, nrep, add_NA = TRUE)








