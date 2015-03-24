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
  z <- array(dim = c(nsite, nspec, nyear))
  lpsi <- array(dim = c(nsite, nspec, nyear))
  psi <- array(dim = c(nsite, nspec, nyear))
  
  for(t in 1:nyear) {
    for (j in 1:nsite) {
      for (i in 1:nspec) {
        
        if (t == 1) { #lpsi = logit of psi
          # just add together colonization
          lgam <- (1 - z0[j, i])* (gam[i])
          # just add together persistence  
          lphi <-  z0[j, i] * (phi[i])
          # put both of them together in the lpsi matrix  
          lpsi[j, i, 1] <- lgam + lphi
          # expit of logit of psi
          psi[j, i, 1] <- expit(lpsi[j, i, 1])
          z[j, i, 1] <- rbinom(1, 1, psi[j, i, 1])
        } else {
          lgam <- (1 - z[j, i, t-1])* (gam[i])
          lphi <- z[j, i, t - 1] * (phi[i])
          lpsi[j, i, t] <- lgam + lphi
          psi[j, i, t] <- expit(lpsi[j, i, t])
          z[j, i, t] <- rbinom(1, 1, psi[j, i, t])
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
#     sim_ymat    sim_ymat    sim_ymat    sim_ymat    sim_ymat    sim_ymat
###----------------------------------------------------------------------------

sim_ymat <- function(sim_list = NULL, z = z,
                     add_NA = TRUE){
  with(sim_list, {
  #makes a nsite, by species, by year, by number of repititions matrix
  x <- array(dim = c(nsite, nspec, nyear, nrep))
  for (j in 1:nsite) {
    for (i in 1:nspec) {
      for (t in 1:nyear) {
        for (k in 1:nrep) {
          x[j, i, t, k] <- rbinom(1, 1, p[i] * z[j, i, t])
        }
      }
    }
  }
  ymat <- array(dim = c(nsite, nspec, nyear))
  for(j in 1:nsite){
    for(i in 1:nspec){
      for(t in 1:nyear){
        ymat[j,i,t] <- sum(x[j,i,t,])
      }
    }
  }
  if(add_NA == TRUE){
    
    years <- floor(runif(5, 1, nyear))
    years <- unique(years)
    for(i in 1:length(years)){
      sites <- floor(runif(15, 1, 100))
      ymat[sites,,years[i]] <- NA
    }
  }
  return(ymat)
  })
}

###############################################################################
###############################################################################
###############################################################################


###----------------------------------------------------------------------------
#     sim_jmat    sim_jmat    sim_jmat    sim_jmat    sim_jmat    sim_jmat
###----------------------------------------------------------------------------

sim_jmat <- function(sim_list = NULL, ymat = ymat){
  with(sim_list, {
  jmat <- array(0, dim = c(nspec, nsite, nyear))
  for( i in 1:nspec){
    for(j in 1:nsite){
      for(t in 1:nyear){
        if(is.na(ymat[j, i, t])==TRUE){
          jmat[i, j, t] <- NA
        }else{
        jmat[i, j, t] <- floor(runif(1,ymat[j, i, t], nrep ))
        }
      }
    }
  }
  return(jmat)
  })
}

###############################################################################
###############################################################################
###############################################################################


j_first <- function(sim_list = NULL){
  with(sim_list, {
    jmat <- array(0, dim = c(nspec, nsite, nyear))
    for( i in 1:nspec){
      for(j in 1:nsite){
        for(t in 1:nyear){
          jmat[i,j,t] <- floor(rnorm(1, 18, 2.5))
        }}}
    jmat[jmat>28] <- 28
    return(jmat)
    
  })
}

y_second <- function(sim_list = NULL, jmat = jmat, z = z){
  ymat <- array(dim = c(nspec, nsite, nyear))
  with(sim_list, {
    for(i in 1:nspec){
      for(j in 1:nsite){
        for(t in 1:nyear){
          ymat[i,j,t] <- rbinom(1, jmat[i, j, t], p[i] * z[j, i, t])
        }
      }
    }
    
    return(ymat)
    
  })
}
