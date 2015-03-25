
###----------------------------------------------------------------------------
#     utility functions   utility functionsutility functions utility functions
###----------------------------------------------------------------------------


logit <- function(x) { 
  log(x/(1 - x))
}

expit <- function(x) {
  exp(x)/(1 + exp(x))
}


###############################################################################
###############################################################################
###############################################################################




###----------------------------------------------------------------------------
#     gen_sim_list      gen_sim_list      gen_sim_list    gen_sim_list 
###----------------------------------------------------------------------------


gen_sim_list <- function(nsite = NULL, nspec = NULL, nyear = NULL,
                         nrep = NULL){
  
  # simulate gamma
  actual_gam <- rbeta(1, 1, 1)
  mu_gam <- logit(actual_gam)
  sd_gam <- rgamma(1, 1)
  log_spec_gam <- rnorm(nspec, mu_gam, sd_gam)
  gam <- expit(log_spec_gam)
  
  # simulate phi
  actual_phi <- rbeta(1, 1, 1)
  mu_phi <- logit(actual_phi)
  sd_phi <- rgamma(1,1)
  log_spec_phi <- rnorm(nspec, mu_phi, sd_phi)
  phi <- expit(log_spec_phi)
  
  # simulate p
  
  actual_p <- runif(1, .1, .3) # keeping it low, because camera traps
  mu_p <- logit(actual_p)
  sd_p <- rgamma(1,1)
  log_spec_p <- rnorm(nspec, mu_p, sd_p)
  p <- expit(log_spec_p)
  
  sim_list <- list(nsite = nsite,
                   nspec = nspec,
                   nyear = nyear,
                   nrep = nrep,
                   p = p,
                   gam = gam,
                   phi = phi,
                   hyperp_sd = list(gam = sd_gam,
                                 phi = sd_phi,
                                 p = sd_p),
                   hyperp_mean = list(gam = actual_gam,
                                 phi = actual_phi,
                                 p = actual_p)
                   )
  
  return(sim_list)
}



###############################################################################
###############################################################################
###############################################################################


###----------------------------------------------------------------------------
#     sim_z   sim_z   sim_z   sim_z   sim_z   sim_z   sim_z   sim_z   sim_z  
###----------------------------------------------------------------------------


sim_z <- function(sim_list = NULL){
  with(sim_list,{
  # initial occupancy states
  # makes a nsite * species matrix for the first season z0
  phi0 <- runif(nspec, 0, 1)
  z0 <- array(dim = c(nsite, nspec))
  for (i in 1:nspec) {
    z0[, i] <- rbinom(nsite, 1, phi0[i])
  }
  # subsequent occupancy
  # makes a nsite by species by year for the following years
  z <- array(dim = c(nspec, nsite, nyear))
  lpsi <- array(dim = c(nspec, nsite, nyear))
  psi <- array(dim = c(nspec, nsite, nyear))
  
  for(t in 1:nyear) {
    for (k in 1:nsite) {
      for (i in 1:nspec) {
        
        if (t == 1) { #lpsi = logit of psi
          # just add together colonization
          lgam <- (1 - z0[k, i])* (logit(gam[i])) #not expit gamma!
          # just add together persistence  
          lphi <-  z0[k, i] * (logit(phi[i]))
          # put both of them together in the lpsi matrix  
          lpsi[i, k, 1] <- lgam + lphi
          # expit of logit of psi
          psi[i, k, 1] <- expit(lpsi[i, k, 1])
          z[i, k, 1] <- rbinom(1, 1, psi[i, k, 1])
        } else {
          lgam <- (1 - z[i, k, t-1])* (gam[i])
          lphi <- z[i, k, t - 1] * (phi[i])
          lpsi[i, k, t] <- lgam + lphi
          psi[i, k, t] <- expit(lpsi[i, k, t])
          z[i, k, t] <- rbinom(1, 1, psi[i, k, t])
        }
      }
    }
  }

  return(z) 
  })
  
}


###############################################################################
###############################################################################
###############################################################################




###----------------------------------------------------------------------------
#     sim_jmat    sim_jmat    sim_jmat    sim_jmat    sim_jmat    sim_jmat
###----------------------------------------------------------------------------

sim_jmat <- function(sim_list = NULL, add_NA = TRUE){
  with(sim_list, {
    jmat <- array(0, dim = c(nspec, nsite, nyear))
      for(k in 1:nsite){
        for(t in 1:nyear){
          jmat[,k,t] <- floor(rnorm(1, 18, 2.5))
        }}
    jmat[jmat>nrep] <- nrep # or less than 0
    jmat[jmat<0] <- 0
    
    if(add_NA == TRUE){
      years <- unique(ceiling(runif(5, 1, nyear)))
      for(i in 1:length(years)){
        sites <- floor(runif(15, 1, nsite))
        jmat[,sites,years[i]] <- NA
      }
    }
    return(jmat)
    
  })
}


###############################################################################
###############################################################################
###############################################################################


###----------------------------------------------------------------------------
#     sim_ymat    sim_ymat    sim_ymat    sim_ymat    sim_ymat    sim_ymat
###----------------------------------------------------------------------------

sim_ymat <- function(sim_list = NULL, jmat = jmat, z = z){
  ymat <- array(dim = c(nspec, nsite, nyear))
  with(sim_list, {
    for(i in 1:nspec){
      for(k in 1:nsite){
        for(t in 1:nyear){
          if(is.na(jmat[i,k,t])==TRUE){
            ymat[i,k,t] <- NA
          }else{
          ymat[i,k,t] <- rbinom(1, jmat[i, k, t], p[i] * z[i, k, t])
          } # end else
        } # end t
      } # end j
    } # end i
    
    return(ymat)
    
  })
}

###############################################################################
###############################################################################
###############################################################################



###----------------------------------------------------------------------------
#     make_zinit      make_zinit    make_zinit    make_zinit    make_zinit
###----------------------------------------------------------------------------

make_zinit <- function(ymat = NULL, simlist = NULL){
  zinit <- array(dim = c(nspec, nsite, nyear))
  for(k in 1:nsite){
    for(i in 1:nspec){
      for(t in 1:nyear){
        zinit[i,k,t] <- ymat[i, k, t]
      }
    }
  }
  zinit[zinit>0] <- 1
  return(zinit)

}

###############################################################################
###############################################################################
###############################################################################

###----------------------------------------------------------------------------
#     sim_matrices sim_matrices sim_matrices sim_matrices sim_matrices 
###----------------------------------------------------------------------------


sim_matrices <- function(sim_list, add_NA = TRUE){
  z <- sim_z(sim_list)
  jmat <- sim_jmat(sim_list, add_NA = add_NA)
  ymat <- sim_ymat(sim_list, jmat, z)
  zinit <- make_zinit(ymat, sim_list)
  
  the_mats <- list(z = z,
                   jmat = jmat,
                   ymat = ymat,
                   zinit = zinit)
  return(the_mats)
}

###############################################################################
###############################################################################
###############################################################################

###----------------------------------------------------------------------------
#     sim_all   sim_all   sim_all   sim_all   sim_all   sim_all
###----------------------------------------------------------------------------

sim_all <- function(nsite = NULL, nspec = NULL, nyear = NULL, nrep = NULL,
                    add_NA = TRUE){
  sim_list <- gen_sim_list(nsite = nsite, nspec = nspec,
                          nyear = nyear, nrep = nrep)
  mats <- sim_matrices(sim_list = sim_list, add_NA = add_NA)
  
  mats$jmat[which(is.na(mats$jmat)==TRUE)] <- 0
  
  data_list <- list(y = mats$ymat, nsite = nsite, nyear = nyear,
                    nspec = nspec, jmat = mats$jmat)
  
  

  
  sims <- list(sim_list = sim_list,
               mats = mats,
               data_list = data_list)
  
  
  return(sims)
}

###############################################################################
###############################################################################
###############################################################################


grab_msd <- function(data = NULL){
  stat <- data$statistics
  ans <- array(0, dim = c(nrow(stat), 2))
  for(i in 1:nrow(stat)){
    ans[i,] <- stat[i, 1:2]
  }
  rownames(ans) <- rownames(stat)
  colnames(ans) <- colnames(stat[,1:2])
  return(ans)
}
grab_quant <- function(data = NULL){
  quant <- data$quantiles
  ans <- array(0, dim = c(nrow(quant), 2))
  for(i in 1:nrow(quant)){
    ans[i,] <- quant[i, c(1,5)]
  }
  rownames(ans) <- rownames(quant)
  colnames(ans) <- colnames(quant)[c(1,5)]
  return(ans)
}

pull_summary <- function(data){
  step_one <- grab_msd(data)
  step_two <- grab_quant(data)
  return(signif(cbind(step_one, step_two), 3))
}
