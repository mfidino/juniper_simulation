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
    
    # p hyperparameters
    p_p ~ dbeta(1, 1)
    mu_p <- log(p_p / (1 - p_p))
    sigma_p ~ dunif(0,10)
    tau_p <- pow(sigma_p, -2)


    

    
    
    
    #### occupancy model, specify the formulation of random effects here
    # species specific random effects
    for (i in 1:(nspec)) {
    gam[i] ~ dnorm(mu_gam, tau_gam)
    phi[i] ~ dnorm(mu_phi, tau_phi)
    lp[i] ~ dnorm(mu_p, tau_p)
    dist2cent_gam[i] ~ dunif(-10, 10)
    dist2cent_phi[i] ~ dunif(-10, 10)
    area_gam[i] ~ dunif(-10, 10)
    area_phi[i] ~ dunif(-10, 10)
    }

    
    
    
    

    for(i in 1:(nspec)){
    for(k in 1:(nsite)){
    z[i, k, 1] ~ dbern(psi_in[i])
    for(t in 2:(nyear)){
    logit(psi[i, k, t]) <-  gam[i] * (1 - z[i, k, t-1]) + 
      (dist2cent_gam[i] * cov[k,1] * (1 - z[i, k, t-1])) +
      (area_gam[i] * cov[k,2] * (1 - z[i, k, t-1])) +
      phi[i] * z[i, k, t-1] + 
      (dist2cent_phi[i] * cov[k,1] *z[i, k, t-1]) +
      (area_phi[i] * cov[k,2] * z[i, k, t-1])
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

   # for(k in 1:(nsite)){
    #  for(t in 1:(nyear)){
   #   sp_rich[k,t] <- sum(z[,k,t])
  #}
#}
  #for(i in 1:(nspec)){
  #  for(p in 1:100){
#beta_est[i,p] <- exp(gam[i] + beta[i] * pred_cov[p]) / (1 + exp(gam[i] + beta[i] * pred_cov[p]))
    #}
  #}



    }