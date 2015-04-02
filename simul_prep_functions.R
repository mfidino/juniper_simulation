
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
                        p_sd = NULL, add_NA = NULL, percent_to_NA = NULL){
  if(missing(nspec)|missing(nsite)|missing(nyear)|missing(nrep)|
     missing(gam)|missing(phi)|missing(p)|missing(gam_sd)|missing(phi_sd)|
     missing(p_sd)|missing(add_NA)|missing(percent_to_NA)){
    stop(c("\n\nSupply all arguments to this function (listed below).\n\n",
           args(make_sim_df)))
  }
  sim_df <- expand.grid(nspec, nsite, nyear, nrep, gam, phi, p, gam_sd,
                        phi_sd, p_sd, add_NA, percent_to_NA)
  colnames(sim_df) <- c("nspec", "nsite", "nyear", "nrep", "gam",
                         "phi", "p", "gam_sd", "phi_sd", "p_sd", "add_NA",
                        "percent_to_NA")
  return(sim_df)
}

###############################################################################
###############################################################################
###############################################################################


###----------------------------------------------------------------------------
#     simulate_from_sim_df    simulate_from_sim_df    simulate_from_sim_df 
###----------------------------------------------------------------------------


simulate_from_sd_df <- function(sim_df = NULL, test = FALSE){
  
  if(test) sim_df <- sim_df[sample(1:nrow(sim_df), 2),]
    
  big_sim <- vector(mode = "list", length = nrow(sim_df))
  
  for( i in 1:nrow(sim_df)){
    
    big_sim[[i]] <- sim_all(nsite = sim_df$nsite[i], 
                            nspec = sim_df$nspec[i],
                            nyear = sim_df$nyear[i],
                            nrep = sim_df$nrep[i],
                            actual_gam = sim_df$gam[i],
                            sd_gam = sim_df$gam_sd[i],
                            actual_phi = sim_df$phi[i],
                            sd_phi = sim_df$phi_sd[i],
                            actual_p = sim_df$p[i], sd_p = sim_df$p_sd[i],
                            add_NA = sim_df$add_NA[i], 
                            percent_to_NA = sim_df$percent_to_NA[i] )
    
    
  }
  
  return(big_sim)
}


###############################################################################
###############################################################################
###############################################################################



###----------------------------------------------------------------------------
#     gen_sim_list      gen_sim_list      gen_sim_list    gen_sim_list 
###----------------------------------------------------------------------------


gen_sim_list <- function(nsite = NULL, nspec = NULL, nyear = NULL,
                         nrep = NULL, actual_gam = NULL, sd_gam = NULL,
                         actual_phi = NULL, sd_phi = NULL,
                         actual_p = NULL, sd_p = NULL,
                         add_NA = NULL, percent_to_NA = NULL){
  
  # simulate gamma
  if(missing(actual_gam)) actual_gam <- rbeta(1, 1, 1)
  if(missing(sd_gam)) sd_gam <- runif(1, 1, 3)
  
  mu_gam <- logit(actual_gam)
  log_spec_gam <- rnorm(nspec, mu_gam, sd_gam)
  gam <- expit(log_spec_gam)
  
  # simulate phi
  if(missing(actual_phi)) actual_phi <- rbeta(1, 1, 1)
  if(missing(sd_phi)) sd_phi <- runif(1, 1, 3)
  
  mu_phi <- logit(actual_phi)
  log_spec_phi <- rnorm(nspec, mu_phi, sd_phi)
  phi <- expit(log_spec_phi)
  
  # simulate p
  
  if(missing(actual_p)) actual_p <- runif(1, .1, .3) # keeping it low, because camera traps
  if(missing(sd_p)) sd_p <- runif(1, 0, 3)
  
  mu_p <- logit(actual_p)
  log_spec_p <- rnorm(nspec, mu_p, sd_p)
  p <- expit(log_spec_p)
  
   # percent to NA
  
  if(missing(add_NA)) add_NA = FALSE
  
  if(missing(percent_to_NA)) percent_to_NA <- 0.2
  
  sim_list <- list(nsite = nsite,
                   nspec = nspec,
                   nyear = nyear,
                   nrep = nrep,
                   gam = gam,
                   phi = phi,
                   p = p,
                   hyperp_sd = list(gam = sd_gam,
                                 phi = sd_phi,
                                 p = sd_p),
                   hyperp_mean = list(gam = actual_gam,
                                 phi = actual_phi,
                                 p = actual_p),
                   add_NA = add_NA,
                   percent_to_NA = percent_to_NA
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
  phi0 <- runif(nspec, .1, .9)
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

sim_jmat <- function(sim_list = NULL){
  with(sim_list, {
    jmat <- array(0, dim = c(nspec, nsite, nyear))
      for(k in 1:nsite){
        for(t in 1:nyear){
          jmat[,k,t] <- floor(rnorm(1, 18, 3))
        }}
    jmat[jmat>nrep] <- nrep # or less than 0
    jmat[jmat<0] <- 0
    
    if(add_NA == TRUE){
      years <- 1:nyear
      for(i in 1:length(years)){
        percent_NA <- ceiling(nsite*percent_to_NA)
        sites <- sample(1:nsite, percent_NA)
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
  
  with(sim_list, {
    ymat <- array(dim = c(nspec, nsite, nyear))
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

make_zinit <- function(ymat = NULL){
  zinit <- ymat

  zinit[zinit>0] <- 1
  return(zinit)

}

###############################################################################
###############################################################################
###############################################################################

###----------------------------------------------------------------------------
#     sim_matrices sim_matrices sim_matrices sim_matrices sim_matrices 
###----------------------------------------------------------------------------


sim_matrices <- function(sim_list){
  z <- sim_z(sim_list)
  jmat <- sim_jmat(sim_list)
  ymat <- sim_ymat(sim_list, jmat, z)
  zinit <- make_zinit(ymat)
  
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
                    actual_gam = NULL, sd_gam = NULL,
                    actual_phi = NULL, sd_phi = NULL,
                    actual_p = NULL, sd_p = NULL, add_NA = NULL,
                    percent_to_NA = NULL){
  sim_list <- gen_sim_list(nsite = nsite, nspec = nspec,
                          nyear = nyear, nrep = nrep, actual_gam = actual_gam,
                          sd_gam = sd_gam, actual_phi = actual_phi,
                          sd_phi = sd_phi, actual_p = actual_p,
                          sd_p = sd_p, add_NA = add_NA,
                          percent_to_NA = percent_to_NA)
  
  mats <- sim_matrices(sim_list = sim_list)
  
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



###----------------------------------------------------------------------------
#     base_file_name    base_file_name    base_file_name    base_file_name
###----------------------------------------------------------------------------



base_file_name <- function(one_from_all_sim = NULL){
  with(one_from_all_sim$sim_list,{
    if(add_NA){
      base_name <- paste("gam(", hyperp_mean$gam,",", hyperp_sd$gam,
                         ")_phi(", hyperp_mean$phi, ",", hyperp_sd$phi,
                         ")_p(", hyperp_mean$p, ",", hyperp_sd$p, 
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

write_diagnostics <- function(mod_mcmc, iter = i, basic_name = basic_name){
  my_line <- gelman.diag(mod_mcmc)$psrf[,2]
  
  if(iter==1){
    write(c( "model", names(my_line)),"gelman_diag_table.txt", sep="\t",
          ncolumns = (length(my_line)+1))
  }
  write(c(basic_name, my_line), "gelman_diag_table.txt", sep = "\t",
        ncolumns = (length(my_line)+1), append=TRUE)
}

###############################################################################
###############################################################################
###############################################################################

###----------------------------------------------------------------------------
#     write_known   write_known   write_known   write_known   write_known   
###----------------------------------------------------------------------------


write_known <- function(one_from_all_sim = all_sim[[i]], iter = i,
                        basic_name = basic_name){
  
  my_line <- c(basic_name, with(one_from_all_sim$sim_list, {
    c(nspec, nsite, nyear, nrep,
      hyperp_mean$gam, hyperp_mean$phi, hyperp_mean$p,
      hyperp_sd$gam, hyperp_sd$phi, hyperp_sd$p,
      add_NA, percent_to_NA)
  }
  ))
  
  column_names <- c("model", "nspec", "nsite", "nyear", "nrep",
                    "mean_gam", "mean_phi", "mean_p",
                    "gam_sd", "phi_sd", "p_sd",
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


write_summary <- function(mod_mcmc = mod_mcmc, iter = i, basic_name = basic_name){
  
  my_line <- pull_summary(summary(mod_mcmc))
  
  column_names <- c("model", names(my_line))
  
  if(iter==1){ # write column names
    write(column_names ,"mcmc_summary.txt", sep="\t",
          ncolumns = length(column_names))
  }
  
  write.table(cbind(basic_name, my_line), "mcmc_summary.txt", append = TRUE,
              sep = "\t", col.names = FALSE, row.names = FALSE)
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

grab_quant <- function(data = NULL){
  quant <- data$quantiles
  ans <- array(0, dim = c(nrow(quant), 3))
  ans <- quant[, c(1,3,5)]
  ans <- signif(ans, 3)
  rownames(ans) <- rownames(quant)
  colnames(ans) <- c("lci", "mode", "uci")
  return(ans)
}



###############################################################################
###############################################################################
###############################################################################




###----------------------------------------------------------------------------
#     pull_summary    pull_summary    pull_summary    pull_summary
###----------------------------------------------------------------------------

pull_summary <- function(data){
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

df_2_array <- function(my_df = NULL){
  require(reshape2)
  my_array <- acast(my_df, species~site~season, value.var = "count")
  dimnames(my_array) <- NULL
  return(my_array)
}


###############################################################################
###############################################################################
###############################################################################




