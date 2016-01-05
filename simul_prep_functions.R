
# functions to be used for simulation analysis

# structure

# 1 = utitlity functions for simulation
# 2 = mcmc summary functions
# 3 = convienience functions




#
#####
######
#1111111111111111111111111111111111111111111111111111111111111111111111111111111
#     utility functions   utility functionsutility functions utility functions
#1111111111111111111111111111111111111111111111111111111111111111111111111111111
######
#####
#

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
#     make_sim_df   make_sim_df   make_sim_df   make_sim_df   make_sim_df 
###----------------------------------------------------------------------------


make_sim_df <- function(nspec = NULL, nsite = NULL, nyear = NULL,
                        nrep = NULL, gam = NULL, phi = NULL,
                        p = NULL, gam_sd = NULL, phi_sd = NULL,
                        p_sd = NULL, beta = NULL, beta_sd = NULL, 
                        beta2 = NULL, beta_sd2 = NULL, inxs = NULL, inxs_sd = NULL,
                        add_NA = NULL, percent_to_NA = NULL,
                        cut_off = NULL, row_replicate = 1){
  # make sure all "really necessary" parameters are in the function
  if(missing(nspec)|missing(nsite)|missing(nyear)|missing(nrep)|
     missing(gam)|missing(phi)|missing(p)|missing(gam_sd)|missing(phi_sd)|
     missing(p_sd)|missing(add_NA)|missing(percent_to_NA)|missing(cut_off)){
    stop(c("\n\nSupply all arguments to this function (listed below).\n\n",
           formals(make_sim_df)))
  }
  # things we may want to change or exclude beta coefs, and inxs
  if(missing(beta)) beta <- beta_sd <- 0
  if(missing(beta2)) beta2 <- beta_sd2 <- 0
  if(missing(inxs)) inxs <- inxs_sd <- 0
  # make a df with all unique combinations of these different parameters
  sim_df <- expand.grid(nspec, nsite, nyear, nrep, gam, phi, p, gam_sd,
                        phi_sd, p_sd, beta, beta_sd, beta2, beta_sd2, inxs, inxs_sd, add_NA, 
                        percent_to_NA, cut_off)
  # cool way to name everything appropriately
  colnames(sim_df) <- head(names(formals(make_sim_df)), -1)
  # make more replicates for simulation
  if(row_replicate>1) sim_df <- sim_df[rep(row.names(sim_df), row_replicate),]
  
  return(sim_df)
}

###############################################################################
###############################################################################
###############################################################################


###----------------------------------------------------------------------------
#     simulate_from_sim_df    simulate_from_sim_df    simulate_from_sim_df 
###----------------------------------------------------------------------------

# wrapper function that takes output from make_sim_df and simulates the data
# also where we can specify the interactions, taken as a vector
simulate_from_sim_df <- function(sim_df = NULL, test = FALSE,
                                 known_inxs = NULL){
  # test is just to determine if things are simulating correctly,
  # it only does 2 simulations instead of all of them
  if(test) sim_df <- sim_df[sample(1:nrow(sim_df), 2),]
  # a list to hold all of the simulations  
  big_sim <- vector(mode = "list", length = nrow(sim_df))
  
  for( i in 1:nrow(sim_df)){
    # SIMULATE
    big_sim[[i]] <- sim_all(nsite = sim_df$nsite[i], 
                            nspec = sim_df$nspec[i],
                            nyear = sim_df$nyear[i],
                            nrep = sim_df$nrep[i],
                            actual_gam = sim_df$gam[i],
                            sd_gam = sim_df$gam_sd[i],
                            actual_phi = sim_df$phi[i],
                            sd_phi = sim_df$phi_sd[i],
                            actual_p = sim_df$p[i], sd_p = sim_df$p_sd[i],
                            actual_beta = sim_df$beta[i],
                            sd_beta = sim_df$beta_sd[i],
                            actual_beta2 = sim_df$beta2[i],
                            sd_beta2 = sim_df$beta_sd[i],
                            actual_inxs = sim_df$inxs[i],
                            sd_inxs = sim_df$inxs_sd[i],
                            add_NA = sim_df$add_NA[i], 
                            percent_to_NA = sim_df$percent_to_NA[i],
                            cut_off = sim_df$cut_off[i],
                            known_inxs = known_inxs)
    
    
  }
  
  return(big_sim)
}


###############################################################################
###############################################################################
###############################################################################



###----------------------------------------------------------------------------
#     gen_sim_list      gen_sim_list      gen_sim_list    gen_sim_list 
###----------------------------------------------------------------------------

# this is a sub-function within the sim_all function
# it takes all of the parameters and data
# and makes the "sim_list" which contains
# the true parameter estimates
# used to simulate the data, number of sites, species, and whatnot.
# the sim_list is used by pretty much every function
# to simulate the the z matrix, j matrix, etc.
gen_sim_list <- function(nsite = NULL, nspec = NULL, nyear = NULL,
                         nrep = NULL, actual_gam = NULL, sd_gam = NULL,
                         actual_phi = NULL, sd_phi = NULL,
                         actual_p = NULL, sd_p = NULL,
                         actual_beta = NULL, sd_beta = NULL,
                         actual_beta2 = NULL, sd_beta2 = NULL,
                         actual_inxs = NULL, sd_inxs = NULL, add_NA = NULL, 
                         percent_to_NA = NULL, cut_off = NULL){
  
  # simulate gamma
  if(missing(actual_gam)) actual_gam <- rbeta(1, 1, 1)
  if(missing(sd_gam)) sd_gam <- runif(1, 1, 3)
  # logit transform this
  mu_gam <- logit(actual_gam)
  # make species specific gammas on logit scale
  log_spec_gam <- rnorm(nspec, mu_gam, sd_gam)
  # convert back to probability scale
  gam <- expit(log_spec_gam)
  
  # simulate phi, if not provided make it up
  if(missing(actual_phi)) actual_phi <- rbeta(1, 1, 1)
  if(missing(sd_phi)) sd_phi <- runif(1, 1, 3)
  # logit transform phi
  mu_phi <- logit(actual_phi)
  # make speices specific phi
  log_spec_phi <- rnorm(nspec, mu_phi, sd_phi)
  # back to probability scale
  phi <- expit(log_spec_phi)
  
  # simulate p, same as the last two
  
  if(missing(actual_p)) actual_p <- runif(1, .1, .3) # keeping it low, because camera traps
  if(missing(sd_p)) sd_p <- runif(1, 0, 3)
  # logit
  mu_p <- logit(actual_p)
  # species specific
  log_spec_p <- rnorm(nspec, mu_p, sd_p)
  # probability
  p <- expit(log_spec_p)
  
  # beta, if missing make it zero!
  if(missing(actual_beta)){
  actual_beta <- 0
  mu_beta <- 0
  sd_beta <- 0
  }
  
  # if we do have the beta
  if(actual_beta>0){
    # logit transform
    mu_beta <- logit(actual_beta)
    # species specific
    log_spec_beta <- rnorm(nspec, mu_beta, sd_beta)
    # back to probability
    beta <- expit(log_spec_beta)
    # make up the scaled covariate we want to use
    covar <- rnorm(nsite)
  }else{mu_beta <- 0 ; beta <- rep(0, nspec); log_spec_beta <- rep(0, nspec); covar <- rep(0, nsite)} 
  # ^^^make everything zero if we don't use it so the later stuff doesnt brick^^^
  
  # beta, this is the same as the last one, look at those comments.
  if(missing(actual_beta2)){
    actual_beta2 <- 0
    mu_beta2 <- 0
    sd_beta2 <- 0
  }
  
  
  if(actual_beta2>0){
    mu_beta2 <- logit(actual_beta2)
    log_spec_beta2 <- rnorm(nspec, mu_beta2, sd_beta2)
    beta2 <- expit(log_spec_beta2)
  }else{mu_beta2 <- 0 ; beta2 <- rep(0, nspec); log_spec_beta2 <- rep(0, nspec)} 
  # no inxs? thats okay, make sd_inxs zero too then.
  if(missing(actual_inxs)){
    sd_inxs <- 0
  }
  # but if you the inxs and their sd
  if(is.numeric(actual_inxs) & is.numeric(sd_inxs)){
    mu_inxs <- 0 # assume mean of zero
    n_inxs <- (nspec^2) - nspec # number of off-diagonal elements in inxs matrix
    log_spec_inxs <- rnorm(n_inxs, mu_inxs, sd_inxs) # simulate the inxs
    inxs <- expit(log_spec_inxs) # back to probability scale
  }else{mu_inxs <- sd_inxs <- 0; n_inxs <- (nspec^2) - nspec; inxs <- rep(0, n_inxs); log_spec_inxs <- rep(0, n_inxs)}
  # ^^^ if no inxs and sd of inxs do this so the later stuff doesnt choke ^^^

  
   # percent to NA
   # if no add_NA, then dont add NA values
  if(missing(add_NA)) add_NA <- FALSE
  # didn't tell me how much you want to NA, then do 20%
  if(missing(percent_to_NA)) percent_to_NA <- 0.2
  
  # cutoff, I can't remember what this does :)
  # I think it means we don't analyze sites at or under cut_off
  # that did not have enough days sampled
  if(missing(cut_off)) cut_off <- 0
  # put our species specific phi values on the diagonal of the phi inxs matrix
  l_phi_mat <- diag(log_spec_phi)
  # put the inxs in, fills column-wise
  l_phi_mat[which(l_phi_mat==0)] <- log_spec_inxs
  # then I guess we transpose it
  # so that it fills row wise instead.
  l_phi_mat <- t(l_phi_mat)
  
  # make the sim_list to return
  # has stuff on logit and probability scale
  sim_list <- list(nsite = nsite,
                   nspec = nspec,
                   nyear = nyear,
                   nrep = nrep,
                   gam = gam,
                   phi = phi,
                   p = p,
                   beta = beta,
                   beta2 = beta,
                   inxs = inxs,
                   hyperp_sd = list(gam = sd_gam, 
                                 phi = sd_phi,
                                 p = sd_p,
                                 beta = sd_beta,
                                 beta2 = sd_beta2,
                                 inxs = sd_inxs),
                   hyperp_mean = list(gam = actual_gam,
                                 phi = actual_phi,
                                 p = actual_p,
                                 beta = actual_beta,
                                 beta2 = actual_beta2,
                                 inxs = actual_inxs),
                   add_NA = add_NA,
                   percent_to_NA = percent_to_NA,
                   l_gam = log_spec_gam,
                   l_phi = log_spec_phi,
                   l_p = log_spec_p,
                   l_beta = log_spec_beta,
                   l_beta2 = log_spec_beta2,
                   l_inxs = log_spec_inxs,
                   cut_off = cut_off,
                   covar = covar,
                   l_phi_mat = l_phi_mat
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
  z <- with(sim_list,{
  # initial occupancy states
  # makes a nsite * species matrix for the first season z0
    
  z <- array(dim = c(nspec, nsite, nyear)) # make z array, species, by site, by season
  logit_psi <- array(dim = c(nspec, nsite, nyear)) # logit_psi array
  psi <- array(dim = c(nspec, nsite, nyear)) # psi on prob scale
  # just randomly making initial probs for occupancy from beta(2,2) distribution
  # may want to change, may not, 
  psi1 <- rbeta(nspec, 2, 2) 
  for (i in 1:nspec) {
    z[i,,1] <- rbinom(nsite, 1, psi1[i])
  }
# subsequent occupancy
# makes a nsite by species by year for the following years

# add together the linear predictor for each other site and season
# put it into logit_psi, transform to a probability, and then
# do a bernoulli trail to fill in the z matrix
  for(t in 2:nyear) {
    for (k in 1:nsite) {
      for (i in 1:nspec) {
          logit_gam <-  (1 - z[i, k, t - 1]) * l_gam[i] + (1 - z[i, k, t - 1]) * l_beta[i] * covar[k]
          logit_phi <- (z[, k, t - 1] %*% l_phi_mat[i,] + l_beta2[i] * covar[k]) * z[i, k, t - 1]
          logit_psi[i, k, t] <- logit_gam + logit_phi
          psi[i, k, t] <- expit(logit_psi[i, k, t])
          z[i, k, t] <- rbinom(1, 1, psi[i, k, t])
        }
      }
    }
  # Note: we could speed this up by pulling the rbinom outside of the loop
  #       so we only have to call it once. We probs should at some point.
  # something like rbinom(nspec*nsite*nyear, 1, psi), to make the
  # z vector and then we put it into the matrix correctly.

  return(z) 
  })
  
}


###############################################################################
###############################################################################
###############################################################################




###----------------------------------------------------------------------------
#     sim_jmat    sim_jmat    sim_jmat    sim_jmat    sim_jmat    sim_jmat
###----------------------------------------------------------------------------

# makes the number of days a camera trap was active matrix (i.e. the j matrix)
sim_jmat <- function(sim_list = NULL){
  with(sim_list, {
    days_array <- array(0, dim = c(nrep, nsite, nyear))
    # the hard coded probability that a camera was active that day
    # calculated from our own data. 
    # note here that we actually sample for 30 days, but then
    # constrain to 28, which is similar to our own camera trap stuff.
    camera_prob <- c(0.45, 0.72, 0.81, 0.84, 0.83, 0.81, 0.82, 0.81, 0.86, 0.88, 
                     0.88, 0.87, 0.87, 0.89, 0.92, 0.94, 0.93, 0.92, 0.90, 0.89,
                     0.89, 0.89, 0.89, 0.88, 0.88, 0.87, 0.86, 0.84, 0.85, 0.85)

# figure out the number of days each camera trap was active
        n_samp <- rbinom(nsite * nyear, length(camera_prob), camera_prob)
        n_samp[n_samp>28] <- 28 # > 28 change to 28
        # randomly sample which days were active for a particular camera.
        # Right now this is done with equal probability.
        n_list <- lapply(n_samp, function(x) sort(c(1, sample(2:28, x-1))))
        
        to_start <- 0
        # for every time a carmera trap was active
        for(i in 1:length(n_list)){
          # put a 1 on the days active on days_array
          days_array[to_start+n_list[[i]]] <- 1
          to_start <- to_start+ nrep # keep trucking through the array
        }
        # reorganize array to line up with our other arrays
        days_array <- aperm(days_array, c(2,3,1))
        # rep days_array by number of species, for jags 
        jmat_expanded <- array(rep(days_array, each = nspec),
                               dim = c(nspec, nsite, nyear, nrep))
    
    # make the jmat 
    jmat <- array(0, dim = c(nspec, nsite, nyear))
    # sum through to give # of days active
    jmat <- apply(jmat_expanded, c(2,3), rowSums)
    
    # this was set up so we could remove sampling in a week
    # if it did not exceed the cutoff, we probably
    # won't use this ever as we are just throwing out data.
    weeks <- rep(1:4, each = 7)
    weeks_to_sum_by <- array(0, dim = c(nspec, nsite, nyear, 4))
    jmat_comparison <- array(0, dim = c(nspec, nsite, nyear))
    
    for(k in 1:nsite){
      for(t in 1:nyear){
        for(j in 1:4){
          weeks_to_sum_by[,k, t, j] <- ifelse(rowSums(jmat_expanded[,k, t, which(weeks==j)])[1]>0,                                              1, 0)
        }     
      }
    }
    jmat_comparison <- apply(weeks_to_sum_by, c(2, 3), rowSums)
    jmat_comparison[which(jmat<cut_off)] <- 0
    

    # to return
    j_list <- list(jmat = jmat, 
                   jmat_comparison = jmat_comparison, 
                   jmat_expanded = jmat_expanded)
    return(j_list)
    
  })
}


###############################################################################
###############################################################################
###############################################################################


###----------------------------------------------------------------------------
#     sim_ymat    sim_ymat    sim_ymat    sim_ymat    sim_ymat    sim_ymat
###----------------------------------------------------------------------------

# this is us making the observed data set
# we sample from the z mat at probability p for a species
# to determine the number of days it was seen on.
sim_ymat <- function(sim_list = NULL, j_list = j_list, z = z){
  
  with(sim_list, {
    # figure out what days it was seen on
    ymat_expanded <- array(rbinom(n = nspec*nsite*nyear*nrep, size = j_list$jmat_expanded, 
                  prob = rep(p * z, nrep)), dim = c(nspec, nsite, nyear, nrep))
    # sum through those days
    ymat <- array(apply(ymat_expanded, c(1,2,3), sum), dim = c(nspec, nsite, nyear))
    
    weeks <- rep(1:4, each = 7)
    # to look at differences between daily or weekly secondary intervals
    # again, probably wont be used.
    ymat_comparison <- array(dim = c(nspec, nsite, nyear))
    sum_by_week <- array(0, dim = c(nspec, nsite, nyear, 4))

    for(j in 1:4){
      sum_by_week[,,,j] <- apply(ymat_expanded[,,,which(weeks==j)],c(1,2,3), sum)
    }
    sum_by_week[sum_by_week>0] <- 1
    
    ymat_comparison <- apply(sum_by_week, c(1,2,3), sum)
    ymat_comparison[which(j_list$jmat_comparison==0)] <- 0
    
    # to return
    y_list <- list(ymat = ymat, ymat_comparison = ymat_comparison, 
                   ymat_expanded = ymat_expanded)
    
    return(y_list)
    
  })
}

###############################################################################
###############################################################################
###############################################################################



###----------------------------------------------------------------------------
#     make_zinit      make_zinit    make_zinit    make_zinit    make_zinit
###----------------------------------------------------------------------------

# to be provided as starting values to the bayesian model
make_zinit <- function(y_list = NULL){
  zinit <- y_list$ymat
  

  zinit[zinit>0] <- 1
  return(zinit)

}

###############################################################################
###############################################################################
###############################################################################

###----------------------------------------------------------------------------
#     make_inxs   make_inxs   make_inxs   make_inxs   make_inxs   make_inxs
###----------------------------------------------------------------------------


# function to appropriately index the interactions
# in the bayesian model.
make_inxs <- function(sim_list = NULL){
  with(sim_list,{
  ns <- nspec
  col_l <- row_l <- col_u <- row_u <-numeric((ns * ns - ns)/2)
  for(i in 1:(ns-1)){
    # algorithm to fill lower columns
    col_l[min(which(col_l ==0)):(min(which(col_l ==0)) + ns-1-i)] <-  rep(i, ns-i)
    # algorithm to fill lower rows
    row_l[min(which(row_l == 0)):(min(which(row_l==0))+ ns-1-i)] <- (i+1):ns
    # algorithm to fill upper rows
    row_u[min(which(row_u == 0)):(min(which(row_u==0))+ ns-(ns+1)+i)] <- 1:i
    # algorithm to fill upper columns
    col_u[min(which(col_u == 0)):(min(which(col_u==0))+ ns-(ns+1)+i)] <- rep(i + 1, i)
    
  }
  
  row_vec <- c(row_l, row_u)
  col_vec <- c(col_l, col_u)

  return(list(rows_vec = row_vec,
              cols_vec = col_vec))})
  
  
}


###############################################################################
###############################################################################
###############################################################################

###----------------------------------------------------------------------------
#     sim_matrices sim_matrices sim_matrices sim_matrices sim_matrices 
###----------------------------------------------------------------------------

# wrapper function to sim all of the matrices,
# this goes inside the sim_all function.
sim_matrices <- function(sim_list = NULL){
  z <- sim_z(sim_list)
  j_list <- sim_jmat(sim_list)
  y_list <- sim_ymat(sim_list, j_list, z)
  
  
# we add the NA values here, post simulation.
      if(sim_list$add_NA == TRUE){
        years <- 1:sim_list$nyear
        for(i in 1:length(years)){
          percent_NA <- ceiling(sim_list$nsite*sim_list$percent_to_NA)
          sites <- sample(1:sim_list$nsite, percent_NA)
          j_list$jmat[,sites,years[i]] <- NA
          j_list$jmat_comparison[,sites, years[i]] <- NA
          y_list$ymat[,sites, years[i]] <- NA
          y_list$ymat_comparison[,sites, years[i]] <- NA
        }
      }
  zinit <- make_zinit(y_list)
  
  the_mats <- list(z = z,
                   j_list = j_list,
                   y_list = y_list,
                   zinit = zinit)

  return(the_mats)
}

###############################################################################
###############################################################################
###############################################################################

###----------------------------------------------------------------------------
#     sim_all   sim_all   sim_all   sim_all   sim_all   sim_all
###----------------------------------------------------------------------------

# this simulates everything, and goes inside simulate_from_sim_df
sim_all <- function(nsite = NULL, nspec = NULL, nyear = NULL, nrep = NULL,
                    actual_gam = NULL, sd_gam = NULL,
                    actual_phi = NULL, sd_phi = NULL,
                    actual_p = NULL, sd_p = NULL,
                    actual_beta = NULL, sd_beta = NULL,
                    actual_beta2 = NULL, sd_beta2 = NULL,
                    actual_inxs = NULL, sd_inxs = NULL, add_NA = NULL,
                    percent_to_NA = NULL, cut_off = NULL,
                    known_inxs = NULL){
  sim_list <- gen_sim_list(nsite = nsite, nspec = nspec,
                          nyear = nyear, nrep = nrep, actual_gam = actual_gam,
                          sd_gam = sd_gam, actual_phi = actual_phi,
                          sd_phi = sd_phi, actual_p = actual_p,
                          sd_p = sd_p, actual_beta = actual_beta,
                          sd_beta = sd_beta, actual_beta2 = actual_beta2,
                          sd_beta2 = sd_beta2, actual_inxs = actual_inxs,
                          sd_inxs = sd_inxs, add_NA = add_NA,
                          percent_to_NA = percent_to_NA,
                          cut_off = cut_off)
# This takes our known_inxs and puts them where they need to be.
  if(length(known_inxs>1)){
    # throw in the inxs
    sim_list$l_phi_mat <- matrix(known_inxs, ncol = sim_list$nspec, nrow = sim_list$nspec)
    diag(sim_list$l_phi_mat) <- sim_list$l_phi # overwrite the diagonals with phis
  }
# simulate the matrices
  mats <- sim_matrices(sim_list = sim_list)
# the jmats cannot have NAs, we know when we never sampled  
  mats$j_list$jmat[which(is.na(mats$j_list$jmat)==TRUE)] <- 0
  mats$j_list$jmat_comparison[which(is.na(mats$j_list$jmat_comparison)==TRUE)] <- 0
  
  # do this if we have a covariate, it will always be more than 1 unique value.
  # if we do. Note that the data_list is for use in jags.
  if(length(unique(sim_list$covar))>1){
    
    data_list <- c(list(y = mats$y_list$ymat, nsite = nsite, nyear = nyear,
                      nspec = nspec, jmat = mats$j_list$jmat, cov = sim_list$covar,
                      n_inxs = nspec^2, only_inxs = (nspec^2 - nspec)),
                   make_inxs(sim_list))
    
    comparison_list <- c(list(y = mats$y_list$ymat_comparison, nsite = nsite, nyear = nyear,
                            nspec = nspec, jmat = mats$j_list$jmat_comparison,
                            cov = sim_list$covar, n_inxs = nspec^2, only_inxs = (nspec^2 - nspec)),
                         make_inxs(sim_list))
  }else{
    data_list <- c(list(y = mats$y_list$ymat, nsite = nsite, nyear = nyear,
                          nspec = nspec, jmat = mats$j_list$jmat,
                        n_inxs = nspec^2, only_inxs = (nspec^2 - nspec)),
                   make_inxs(sim_list))

  
    comparison_list <- c(list(y = mats$y_list$ymat_comparison, nsite = nsite, nyear = nyear,
                          nspec = nspec, jmat = mats$j_list$jmat_comparison,
                          n_inxs = nspec^2, only_inxs = (nspec^2 - nspec)),
                         make_inxs(sim_list))
  }
  

  
  sims <- list(sim_list = sim_list,
               mats = mats,
               data_list = data_list,
               comparison_list = comparison_list)
  
  
  return(sims)
}

###############################################################################
###############################################################################
###############################################################################



###----------------------------------------------------------------------------
#     base_file_name    base_file_name    base_file_name    base_file_name
###----------------------------------------------------------------------------


# wrapper function to write the file name with
# the simulated values.
base_file_name <- function(one_from_all_sim = NULL){
  with(one_from_all_sim$sim_list,{
    if(add_NA){
      base_name <- paste("gam(", hyperp_mean$gam,"_", hyperp_sd$gam,
                         ")_phi(", hyperp_mean$phi, "_", hyperp_sd$phi,
                         ")_p(", hyperp_mean$p, "_", hyperp_sd$p, 
                         ")_NA(", percent_to_NA, ")",
                         sep = "")
    }else{
      base_name <-paste("gam(", hyperp_mean$gam,",", hyperp_sd$gam,
                        ")_phi(", hyperp_mean$phi, ",", hyperp_sd$phi,
                        ")_p(", hyperp_mean$p, ",", hyperp_sd$p, 
                        ")_NA(", "0.0)",
                        sep = "")
    }
    return(base_name)
  })
  
}


###############################################################################
###############################################################################
###############################################################################


###----------------------------------------------------------------------------
#     write_mcmc_matrix   write_mcmc_matrix   write_mcmc_matrix   
###----------------------------------------------------------------------------

# write the mcmc_matrix. 
write_mcmc_matrix <- function(mod_mcmc = NULL, basic_name = NULL){
  mod_mcmc_path <- paste("C:/simulations/dynamic_occupancy/mcmc_matrix/mcmc_matrix_",
                         basic_name, ".txt", sep = "")
  write.table(as.matrix(mod_mcmc, chains = TRUE), mod_mcmc_path,
              row.names = FALSE, sep = "\t")
  
}

###############################################################################
###############################################################################
###############################################################################


###----------------------------------------------------------------------------
#     write_diagnostics   write_diagnostics   write_diagnostics  
###----------------------------------------------------------------------------

# note, this is a sub-function that works within a for loop
# does a gelman diagnostic each iterations and we
# append that into a table.
write_diagnostics <- function(mod_mcmc, iter = i, basic_name = basic_name){
  my_line <- gelman.diag(mod_mcmc)$psrf[,2]
  
  if(iter==1){
    write(c( "model", names(my_line)),"diagnostic_table.txt", sep="\t",
          ncolumns = (length(my_line)+1))
  }
  write(c(paste(basic_name,"iter", iter, sep = "_"), my_line), "diagnostic_table.txt", sep = "\t",
        ncolumns = (length(my_line)+1), append=TRUE)
}

###############################################################################
###############################################################################
###############################################################################

###----------------------------------------------------------------------------
#     write_known   write_known   write_known   write_known   write_known   
###----------------------------------------------------------------------------

# also used within an iteration, basically writes all of our relevant info
# from the sim_list. Not the z matrix though.
write_known <- function(one_from_all_sim = all_sim[[i]], iter = i,
                        basic_name = basic_name, mod_mcmc = mod_mcmc){
  
  my_line <- c(paste(basic_name,"iter", iter, sep = "_"), with(one_from_all_sim$sim_list, {
    c(nspec, nsite, nyear, nrep,
      gam, logit(hyperp_mean$gam), logit(hyperp_mean$phi),
      p, hyperp_mean$gam, hyperp_mean$p, hyperp_mean$phi,
      phi, rep(NA, nspec), 
      hyperp_sd$gam, hyperp_sd$p, hyperp_sd$phi,
      add_NA, percent_to_NA)
  }
  ))
  
  column_names <- c("model", "nspec", "nsite", "nyear", "nrep",
                    names(gelman.diag(mod_mcmc)$psrf[,2]),
                    "add_NA", "percent_to_NA")
  
  
  if(iter==1){ # write column names
  
    write(column_names ,"known_values.txt", sep="\t",
          ncolumns = length(column_names))
  }
  
  write(my_line, "known_values.txt", append = TRUE,
        sep = "\t", ncolumns = length(column_names))
}


###############################################################################
###############################################################################
###############################################################################



###----------------------------------------------------------------------------
#     write_summary   write_summary   write_summary   write_summary 
###----------------------------------------------------------------------------

# grabs relavant summary info from the posterior distributions,
# again, used within a loop.
write_summary <- function(mod_mcmc = mod_mcmc, iter = i, basic_name = basic_name){
  
  my_line <- pull_summary(summary(mod_mcmc))
  
  column_names <- c("model", names(my_line))
  
  if(iter==1){ # write column names
    write(column_names ,"mcmc_summary.txt", sep="\t",
          ncolumns = length(column_names))
  }
  
  write.table(cbind(paste(basic_name,"iter", iter, sep = "_"), my_line), "mcmc_summary.txt", append = TRUE,
              sep = "\t", col.names = FALSE, row.names = FALSE)
}


###############################################################################
###############################################################################
###############################################################################



###----------------------------------------------------------------------------
#     batch_analyze   batch_analyze   batch_analyze   batch_analyze 
###----------------------------------------------------------------------------

# big loop to batch analyze the models.

batch_analyze <- function(all_sim = NULL, params = NULL,
                          n_chains = NULL, adapt_steps = NULL,
                          burn_in = NULL, sample_steps = NULL,
                          thin_steps = NULL, make_comparisons = NULL,
                          model = NULL, m_name = NULL){
  
  for(i in 1:length(all_sim)){
    
    print(paste("Analyzing", i, "of", length(all_sim), "simulations", sep = " "))
    if(make_comparisons) print("Analyzing data_list")
   
#  # generate initial values
#     inits <- function(chain){
#       gen_list <- function(chain = chain){
#         list( 
#           z = all_sim[[i]]$mats$zinit,
#           gam = rbeta(1, 1, 1),
#           sigma_gam = runif(1, 0, 5),
#           p_phi = rbeta(1, 1, 1),
#           sigma_phi = runif(1, 0, 5),
#           p_p = rbeta(1, 1, 1),
#           .RNG.name = switch(chain,
#                              "1" = "base::Wichmann-Hill",
#                              "2" = "base::Marsaglia-Multicarry",
#                              "3" = "base::Super-Duper",
#                              "4" = "base::Mersenne-Twister",
#                              "5" = "base::Wichmann-Hill",
#                              "6" = "base::Marsaglia-Multicarry",
#                              "7" = "base::Super-Duper",
#                              "8" = "base::Mersenne-Twister"),
#           .RNG.seed = sample(1:1e+06, 1)
#           )
#       }
#       return(switch(chain,           
#         "1" = gen_list(chain),
#         "2" = gen_list(chain),
#         "3" = gen_list(chain),
#         "4" = gen_list(chain),
#         "5" = gen_list(chain),
#         "6" = gen_list(chain),
#         "7" = gen_list(chain),
#         "8" = gen_list(chain)
#         )
#         )
#     }
    inits <- function(chain){
      gen_list <- function(chain = chain){
        list( 
          z = all_sim[[i]]$mats$zinit,
          gam = runif(3, -3, 3),
          spe_phi = runif(3, -3, 3),
          lp = runif(3, -3, 3),
          tau_int = runif(1, 0.01, 5),
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

    
# run the jags model.
    mod_mcmc2 <- as.mcmc.list(run.jags( model= model , 
                                       monitor=params , 
                                       data=all_sim[[i]]$data_list ,  
                                       inits=inits , 
                                       n.chains=n_chains ,
                                       adapt=adapt_steps ,
                                       burnin=burn_in , 
                                       sample=ceiling(sample_steps / n_chains) ,
                                       thin=thin_steps ,
                                       summarise=FALSE ,
                                       plots=FALSE,
                                       method = "parallel"))
    
    mod_mcmc_path <- paste("C:/simulations/dynamic_occupancy/mcmc_matrix/mcmc_matrix_",
                           m_name,i, ".txt", sep = "")
    write.table(as.matrix(mod_mcmc, chains = TRUE), mod_mcmc_path,
                row.names = FALSE, sep = "\t")
saveRDS(all_sim[[i]], paste("C:/simulations/dynamic_occupancy/knowns/all_sim_",
                                i, ".txt", sep = ""))
# currently turned this off, we should probably save this stuff though. 
#     # writing out file stuff
#     basic_name <- base_file_name(all_sim[[i]])
#     
#     write_diagnostics(mod_mcmc, iter = i, basic_name = basic_name)
#     
#     write_summary(mod_mcmc, iter = i, basic_name = basic_name)
#     
#     write_known(one_from_all_sim = all_sim[[i]], iter = i, basic_name = basic_name,
#                 mod_mcmc = mod_mcmc)

    # if you are comparing inxs to non_inxs model, run it again
    # with beta only model.
    if(make_comparisons){
      print("Analyzing with beta_model")
      
      inits_beta <- function(chain){
        gen_list <- function(chain = chain){
          list( 
            z = all_sim[[i]]$mats$zinit,
            gam = runif(3, -3, 3),
            phi = runif(3, -3, 3),
            lp = runif(3, -3, 3),
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
      mod_mcmc_compare <- as.mcmc.list(run.jags( model="beta_model.R" , 
                                         monitor=params , 
                                         data=all_sim[[i]]$data_list,  
                                         inits=inits_beta , 
                                         n.chains=n_chains ,
                                         adapt=adapt_steps ,
                                         burnin=burn_in , 
                                         sample=ceiling(sample_steps / n_chains) ,
                                         thin=thin_steps ,
                                         summarise=FALSE ,
                                         plots=FALSE,
                                         method = "parallel"))
      
      
      mod_mcmc_path <- paste("C:/simulations/dynamic_occupancy/mcmc_matrix/mcmc_beta_matrix_",
                             m_name,i, ".txt", sep = "")
      write.table(as.matrix(mod_mcmc_compare, chains = TRUE), mod_mcmc_path,
                  row.names = FALSE, sep = "\t")

#       compare_i <- 2
#       write_diagnostics(mod_mcmc_compare, iter = compare_i, basic_name = basic_name)
#       
#       write_summary(mod_mcmc_compare, iter = compare_i, basic_name = basic_name)
#       
#       write_known(one_from_all_sim = all_sim[[i]], iter = compare_i, 
#                   basic_name = basic_name,
#                   mod_mcmc = mod_mcmc_compare)
    }
  }
}


###############################################################################
###############################################################################
###############################################################################


#
#####
######
#222222222222222222222222222222222222222222222222222222222222222222222222222222
# MCMC summary functions    MCMC summary functions    MCMC summary functions
#222222222222222222222222222222222222222222222222222222222222222222222222222222
######
#####
#


###----------------------------------------------------------------------------
#     grab_msd    grab_msd    grab_msd    grab_msd    grab_msd    grab_msd
###----------------------------------------------------------------------------

# collects the mean from a the mcmc post summary
grab_msd <- function(data = NULL){
  stat <- data$statistics
  ans <- array(0, dim = c(nrow(stat), 1))
    ans <- stat[, 1]
  
  ans <- signif(ans, 3)
  ans <- data.frame(rownames(stat), ans)
  colnames(ans) <- c("parameter", "mean")
  return(ans)
}



###############################################################################
###############################################################################
###############################################################################


###----------------------------------------------------------------------------
#     grab_quant    grab_quant    grab_quant    grab_quant    grab_quant
###----------------------------------------------------------------------------

# grab upper, lower, and median from post mcmc summary

grab_quant <- function(data = NULL){
  quant <- data$quantiles
  ans <- array(0, dim = c(nrow(quant), 3))
  ans <- quant[, c(1,3,5)]
  ans <- signif(ans, 3)
  rownames(ans) <- rownames(quant)
  colnames(ans) <- c("lci", "median", "uci")
  return(ans)
}



###############################################################################
###############################################################################
###############################################################################




###----------------------------------------------------------------------------
#     pull_summary    pull_summary    pull_summary    pull_summary
###----------------------------------------------------------------------------
# used in an above summary function.
pull_summary <- function(data){
  if(class(data) == "mcmc.list") data <- summary(data)
  step_one <- grab_msd(data)
  step_two <- grab_quant(data)
  return(data.frame(step_one, step_two))
}

###############################################################################
###############################################################################
###############################################################################


#
#####
######
#333333333333333333333333333333333333333333333333333333333333333333333333333333
# Convienience functions    Convienience functions    Convienience functions
#333333333333333333333333333333333333333333333333333333333333333333333333333333
######
#####
#




###----------------------------------------------------------------------------
#     array_2_df    array_2_df    array_2_df    array_2_df    array_2_df
###----------------------------------------------------------------------------

# change and array to a dataframe.
array_2_df <- function(my_array = NULL){
  require(plyr)
  if(length(dim(my_array)) != 3){
    stop("Error: You did not give this function a 3-dimensional array.")
  }
  
  my_df <- adply(my_array, c(1,2,3))
  colnames(my_df) <- c("species", "site", "season", "count")
  return(my_df)
  
}

###############################################################################
###############################################################################
###############################################################################



###----------------------------------------------------------------------------
#     df_2_array    df_2_array    df_2_array    df_2_array    df_2_array
###----------------------------------------------------------------------------

# change the dataframe back to an array
df_2_array <- function(my_df = NULL){
  require(reshape2)
  my_array <- acast(my_df, species~site~season, value.var = "count")
  dimnames(my_array) <- NULL
  return(my_array)
}


###############################################################################
###############################################################################
###############################################################################




