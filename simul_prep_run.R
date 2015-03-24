#########################
#
#
# Working example of dynamic multi=-species occupancy model
#
# 1/26/2014
#
#

########################
# Current functionality#
########################

# This model is set for species level random effects in regards to 
# occupancy, colonization, persistence, and detection probability.

# Additionally, I have included one covariate for both colonization and persistence
# To code for persistence level effect use " * z[j, i, t - 1]" while  colonization is
# by "* 1 - z[j, i, t-1]"

######################################################
# more things to do on the model:
######################################################

# Include more than just random species level effect in detection
# spatially explicit rescue effects
# interaction between species
# nsite specific estimation of species richness / what species are there
# allow for number of revisits / nsite to change depending on how many active days are there



# utility functions

logit <- function(x) { 
    log(x/(1 - x))
}

expit <- function(x) {
    exp(x)/(1 + exp(x))
}


# simulate the data


nsite <- 100
nspec <- 2
nyear <- 5
nrep <- 28



# community level hyperparameters



# current settings 


# gam= colonization
p_gam= 0.3
mugam<- logit(p_gam)
sdgam <- 2

# phi = persistence
p_phi <- 0.6
muphi <- logit(p_phi)
sdphi <- 1

# generate observed data
p_p <- 0.25
mup <- logit(p_p)
sdp <- 1
lp <- rnorm(nspec, mup, sdp)
p <- expit(lp)






# species specific random effects
set.seed(1)  # for reproducibility
gam <- rnorm(nspec, mugam, sdgam)
set.seed(1008)
phi <- rnorm(nspec, muphi, sdphi)

sim_list <- list(nsite = nsite,
                 nspec = nspec,
                 nyear = nyear,
                 nrep = nrep,
                 p = p,
                 gam = gam,
                 phi = phi)


# make z matrix, with NA values

z <- sim_z(sim_list)

 # places to mess with the code, add NA and missing data, etc. 
# Condense to binomial, j matrix
# start with j same
# then change j

ymat <- sim_ymat(sim_list, z, add_NA =TRUE)

# put data togther, we should probably figure out a way to keep all covariate data in one
# spot so we don't have to specify it in the data here one by one.

jmat <- sim_jmat(sim_list, ymat)



win.data <- list(x=x, nrep=nrep, nsite=nsite, nspec = nspec, nyear = nyear, 
                 cov = data.frame(isol_cov))
# allow for temporal component for covariates

# initial values, if we did observed a species at a nsite put z = 1

zinit <- array(dim = c(nsite, nspec, nyear))
for (j in 1:nsite) {
  for (i in 1:nspec) {
    for (t in 1:nyear) {
      zinit[j, i, t] <- max(x[j, i, t, ]) # change to one if greater than.
    }
  }
} #basically get max of x for this

inits <- function() { # specify initial values
  list(p_gam = runif(1, 0, 1), p_phi = runif(1, 0, 1), 
       sigmaphi = runif(1, 0, 1), sigmap = runif(1, 0, 3), 
       sigmagam = runif(1, 0, 3), p_isolp = runif(1, 0 ,1),
       p_isolg = runif(1, 0 ,1), z = zinit)
}

# keep track of the paramters
params <- c( "gam", "phi", "z", "isolg", "isolp")


#write out the model

sink("dynamic_random_area.txt")
cat("
    model{
    #### priors
    
    # gam hyperparameters
    p_gam ~ dbeta(1, 1)
    mugam <- log(p_gam / (1 - p_gam))
    sigmagam ~ dunif(0, 10)
    taugam <- (1 / (sigmagam * sigmagam))
    
    # phi hyperparameters
    p_phi ~ dbeta(1, 1)
    muphi <- log(p_phi / (1 - p_phi))
    sigmaphi~dunif(0,10)
    tauphi<-1/(sigmaphi*sigmaphi)
    
    # p hyperparameters
    p_p ~ dbeta(1, 1)
    mup <- log(p_p / (1 - p_p))
    sigmap ~ dunif(0,10)
    taup <- (1 / (sigmap * sigmap))



    # isolation hyperparamters for gamma
    p_isolg ~ dbeta(1,1)
    muisolg <- log(p_isolg / (1 - p_isolg))
    sigmaisolg ~ dunif(0, 10)
    tauisolg <- 1 / (sigmaisolg*sigmaisolg)

    # isolation hyperparamters for phi
    p_isolp ~ dbeta(1,1)
    muisolp <- log(p_isolp / (1 - p_isolp))
    sigmaisolp ~ dunif(0, 10)
    tauisolp <- 1 / (sigmaisolp*sigmaisolp)
    

    #### occupancy model, specify the formulation of random effects here
    # species specific random effects
    for (i in 1:(nspec)) {
    phi0[i] ~ dbeta(1, 1)
    gam[i] ~ dnorm(mugam, taugam)
    phi[i] ~ dnorm(muphi, tauphi)
    isolp[i] ~ dnorm(muisolp, tauisolp)
    isolg[i] ~ dnorm(muisolg, tauisolg)
}




    
    # occupancy states
    for (j in 1:nsite) {
    for (i in 1:nspec) {
    z0[j, i] ~ dbern(phi0[i])

    
    logit(psi[j, i, 1]) <-gam[i] * (1- z0[j, i]) + (isolg[i] * cov[j, ]) * (1- z0[j, i])
                      + z0[j, i] * phi[i] + z0[j, i] *(isolp[i] * cov[j, ])
    z[j, i, 1] ~ dbern(psi[j, i, 1]) 
    for (t in 2:nyear) {
    logit(psi[j, i, t]) <- (1 - z[j, i, t-1]) * gam[i] + (1 - z[j, i, t-1]) *(isolg[i] * cov[j, ])
                            + z[j, i, t-1] * phi[i]+ z[j, i, t-1] *  (isolp[i] * cov[j, ])
    z[j, i, t] ~ dbern(psi[j, i, t])
    }
    }
    }
    
    
    
    
    #### detection model # this needs to be reworked
    for(i in 1:nspec){ 
    lp[i] ~ dnorm(mup, taup) 
    p[i] <- (exp(lp[i])) / (1 + exp(lp[i])) #This is the spot to fix correct?
    }
    
    
    #### observation model
    for (t in 1:nyear){
    for (j in 1:nsite){
    for (i in 1:nspec){
    
    mu[j, i, t] <- z[j, i, t] * p[i] 
    for (k in 1:nrep){
    x[j, i, t, k] ~ dbern(mu[j, i, t])
    }
    }
    }
    }
    }
    ", fill=TRUE)
sink()


#JUST USE JAGS


ocmod <- jags.model(file = "dynamic_random_area.txt", 
                    inits = inits, data = win.data, n.chains = 2)
nburn <- 2000
update(ocmod, n.iter = nburn)
out <- coda.samples(ocmod, n.iter = 8000, variable.names = params)



nsite <- 100
nspec <- 3
nyear <- 4

test2 <- data.frame(sum_out$quantiles)

to_keep <- grep("z", row.names(test2))

richne <- test2[to_keep,]



yr_one <- matrix(richne[1:300,1], ncol = nspec, nrow = nsite)

for()