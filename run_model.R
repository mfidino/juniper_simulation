inits <- function(){ # Must be = instead of <-
  list( 
        z = sims$mats$zinit,
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



ocmod <- jags.model(file = "intercept_model.txt",
                    inits = inits, data = sims$data_list, n.chains = 2)
nburn <- 10000
update(ocmod, n.iter = nburn)

m1 <- coda.samples(ocmod, n.iter = 10000, variable.names = params)

m1_sum <- summary(m1)

pull_summary(m1_sum)
