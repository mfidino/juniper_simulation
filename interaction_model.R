model{
    #### priors
    
    for(i in 1:nspec){
    psi_in[i] ~ dbeta(1, 1)
    }


# gam hyperparameters
    p_gam ~ dbeta(1, 1)
    mu_gam <- log(p_gam / (1 - p_gam))
    sigma_gam ~ dunif(0, 10)
    tau_gam <- pow(sigma_gam, -2)
    
    # phi hyperparameters
    p_phi ~ dbeta(1, 1)
    mu_phi <- log(p_phi / (1 - p_phi))
    sigma_phi~dunif(0,10)
    tau_phi<- pow(sigma_phi, -2)
    
    # interaction hyperparameters
    p_int ~ dbeta(2,2)
    mu_int <- log(p_int / (1 - p_int))
    sigma_int~dunif(0, 10)
    tau_int <- pow(sigma_int, -2)
    
    # p hyperparameters
    p_p ~ dbeta(1, 1)
    mu_p <- log(p_p / (1 - p_p))
    sigma_p ~ dunif(0,10)
    tau_p <- pow(sigma_p, -2)


    

    
    
    
    #### occupancy model, specify the formulation of random effects here
    # species specific random effects
    
    for (i in 1:(nspec)){
      gam[i] ~ dnorm(mu_gam, tau_gam)
      lp[i] ~ dnorm(mu_p, tau_p)
      beta[i] ~ dunif(-10, 10)
      beta2[i] ~ dunif(-10, 10)
      spe_phi[i] ~ dnorm(mu_phi, tau_phi)
    }
    
    
      for(k in 1:(only_inxs)){
      int_phi[k] ~ dnorm(mu_int, tau_int)
      }
    
      for(i in 1:(n_inxs)){
          phi[rows_vec[i],cols_vec[i]] <- ifelse(equals(rows_vec[i], cols_vec[i]),
                                                 spe_phi[rows_vec[i]], ifelse(equals(int_phi[inxs_vec[i]], 0),
                                                                              0, int_phi[inxs_vec[i]]))
        }
      
      
      


    for(i in 1:(nspec)){
    for(k in 1:(nsite)){
    z[i, k, 1] ~ dbern(psi_in[i])
    for(t in 2:(nyear)){
    logit(psi[i, k, t]) <-  gam[i] * (1 - z[i, k, t-1]) + 
      (beta[i] * cov[k] * (1 - z[i, k, t-1])) +
       (inprod(phi[i,], z[, k, t-1])* z[i, k, t-1]) + (beta2[i] * cov[k] *z[i, k, t-1]) 
    z[i, k, t] ~ dbern(psi[i, k, t])
    }
    }
    }
    
    
    
    
    #### detection model # this needs to be reworked
    for(i in 1:(nspec)){ 
    p[i] <- (exp(lp[i])) / (1 + exp(lp[i])) #This is the spot to fix correct?
    }
    
    
    #### observation model
    for(i in 1:(nspec)){
    for(k in 1:(nsite)){
    for(t in 1:(nyear)){
    mu[i, k, t] <- z[i, k, t] * p[i]
    y[i, k, t] ~ dbin(mu[i, k, t], jmat[i, k, t])
    }
    }
    }





    }