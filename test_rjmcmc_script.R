
data.list <- list(
  nsite = 80,
  nyear = 6,
  nspec = 3,
  jmat = all_sim[[1]]$data_list$jmat,
  cov = all_sim[[1]]$data_list$cov,
  y = all_sim[[1]]$data_list$y
)

inits = list( z = all_sim[[1]]$mats$zinit,
              m = 2)

inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = all_sim[[i]]$mats$zinit,
      m = 2,
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}




m1 <- as.matrix(as.mcmc.list(mod_mcmc))

m1_mean <- apply(as.matrix(as.mcmc.list(mod_mcmc)), 2, mean )
m1_sd <- apply(as.matrix(as.mcmc.list(mod_mcmc)), 2, sd )
m1_log <- apply(as.matrix(as.mcmc.list(mod_mcmc))[,27:29], 2, logit )
apply(m1_log, 2, sd )

cbind(t(t(m1_mean)), t(t(m1_sd)))

m2_mean <- apply(as.matrix(as.mcmc.list(mod_mcmc2)), 2, mean )
m2_sd <- apply(as.matrix(as.mcmc.list(mod_mcmc2)), 2, sd )

m1_log <- apply(as.matrix(as.mcmc.list(mod_mcmc2))[,21:23], 2, logit )
apply(m1_log, 2, sd )
 
cbind(t(t(m2_mean)), t(t(m2_sd)))

params <- c("m", 
            "psi_in1", 
            "psi_in2", 
            "gam1", 
            "gam2",
            "beta1", 
            "beta2",
            "phi1",
            "phi2",
            "lp1",
            "lp2")



test3 <- run.jags( model="test_rjmcmc.R" , 
                       monitor= "m" , 
                       data=data.list ,  
                       inits=inits , 
                       n.chains=7 ,
                       adapt=80 ,
                       burnin=0 , 
                       sample=ceiling(700 / 7) ,
                       thin=thin_steps ,
                       summarise=FALSE ,
                       plots=FALSE,
                       method = "parallel")

t_mean <- apply(as.matrix(as.mcmc.list(test3)), 2, mean)
t_sd <- apply(as.matrix(as.mcmc.list(test3)), 2, sd )
 
 
 
 
test_out <- as.matrix(as.mcmc.list(test3), chain = TRUE)
table(test_out[,2])
 
 
 
 










