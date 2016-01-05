model{
    #### priors
    
    for(i in 1:nspec){
    psi_in[i] ~ dbeta(1, 1)
    }

    
    # interaction hyperparameters
    #p_int ~ dbeta(2,2)
    #mu_int <- log(p_int / (1 - p_int))
    sigma_int <- 1 / sqrt(tau_int)
    tau_int ~ dgamma(1, 1)
    
    
    #### occupancy model, specify the formulation of random effects here
    # species specific random effects
    
    for (i in 1:(nspec)){
      gam[i] ~ dunif(-10, 10)
      lp[i] ~ dunif(-10, 10)
      beta[i] ~ dunif(-10, 10)
      beta2[i] ~ dunif(-10, 10)
      spe_phi[i] ~ dunif(-10, 10)
    }
    
    
      for(k in 1:(only_inxs)){
      int_phi[k] ~ dnorm(0, tau_int)
      }
      
      for(i in 1:nspec){
        phi[i,i] <- spe_phi[i]
      }
      for(i in 1:(n_inxs - nspec)){
        phi[rows_vec[i], cols_vec[i]] <- int_phi[i]
      }
#      for(i in 1:(n_inxs)){
#          phi[rows_vec[i],cols_vec[i]] <- equals(rows_vec[i], cols_vec[i]) * spe_phi[rows_vec[i]]+
#                                          (rows_vec[i] != cols_vec[i]) * int_phi[inxs_vec[i]]
#            
#           # ifelse(equals(rows_vec[i], cols_vec[i]),
#                 #                                spe_phi[rows_vec[i]], int_phi[inxs_vec[i]])
#        }
#      
      
      


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
  var phi[nspec,nspec]