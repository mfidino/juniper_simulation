#########################
#
#
# simul_prep_run for simulation work
#
# 3/24/2015
#
# written by M. Fidino



source("simul_prep_function.R")



# simulate the data


nsite <- 100
nspec <- 2
nyear <- 5
nrep <- 28


sim_list <- gen_sim_list(nsite, nspec, nyear, nrep)


# make z matrix

mats <- sim_matrices(sim_list)








