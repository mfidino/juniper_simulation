inits <- function(){ # Must be = instead of <-
  list( 
        z = big_sim[[10]]$mats$zinit,
         p_gam = rbeta(1,1,1),
        sigma_gam = runif(1, 0, 5),
        p_phi = rbeta(1,1,1),
        sigma_phi = runif(1, 0, 5),
        p_p = rbeta(1, 1, 1)
  )
}


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

n_chains = 2
adapt_steps = 2000
burn_in = 15000
sample_steps = 15000
thin_steps = 1

test3 <- run.jags( model="intercept_model.txt" , 
          monitor=params , 
          data=big_sim[[10]]$data_list ,  
          inits=inits , 
          n.chains=n_chains ,
          adapt=adapt_steps ,
          burnin=burn_in , 
          sample=sample_steps ,
          thin=thin_steps ,
          summarise=FALSE ,
          plots=FALSE )

test3 <- as.mcmc.list(test3)

test_sum3 <- summary(test3)
pull_summary(test_sum3)

caterplot(test3, "gam")
caterpoints(big_sim[[10]]$sim_list$gam)
abline(v = 1)
