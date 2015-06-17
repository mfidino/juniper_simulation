model{


m ~ dcat(p_model[1:2])
p_model[1] <- 0.51
p_model[2] <- 0.5



for(i in 1:(nspec)){
  for(k in 1:(nsite)){
    z[i, k, 1] ~ dbern(psi_in[i])
    for(t in 2:(nyear)){
      logit(psi[i, k, t]) <-  (gam[i] * (1 - z[i, k, t-1]) + 
        (beta[i] * cov[k] * (1 - z[i, k, t-1])) +
        inprod(phi[i, ], z[, k, t-1]) * z[i, k, t-1])
      z[i, k, t] ~ dbern(psi[i, k, t])
    }
  }
}
psi_in <- equals(m, 1)*psi_in1[,1] + equals(m, 2) * psi_in1[,1]




# for first model  
psi_in1[1,1] ~ dnorm(0.805, pow(0.04375, -2))
psi_in1[2,1] ~ dnorm(0.637, pow( 0.0731, -2))
psi_in1[3,1] ~ dnorm(0.709, pow(0.0501, -2))
psi_in1[1,2] ~ dnorm(0.805, pow(0.04375, -2))# pseudo
psi_in1[2,2] ~ dnorm(0.637, pow( 0.0731, -2))# pseudo
psi_in1[3,2] ~ dnorm(0.709, pow(0.0501, -2))# pseudo

#for second model
psi_in2[1,1] ~ dnorm(0.806, pow(0.0435, -2))# psuedo
psi_in2[2,1] ~ dnorm(0.66, pow(0.0762, -2))#psuedo
psi_in2[3,1] ~ dnorm(0.709, pow(0.0499, -2))#psuedo
psi_in2[1,2] ~ dnorm(0.806, pow(0.0435, -2))
psi_in2[2,2] ~ dnorm(0.66, pow(0.0762, -2))
psi_in2[3,2] ~ dnorm(0.709, pow(0.0499, -2))

  
gam <- equals(m,1) * gam1[,1] + equals(m,2) * gam1[,1] 

# first model
gam1[1,1] ~ dnorm(2.12, pow(0.4698, -2))
gam1[2,1] ~ dnorm(-1.89 , pow(1.4825, -2))
gam1[3,1] ~ dnorm(2.52 , pow(0.4313, -2))
gam1[1,2] ~ dnorm(2.12, pow(0.4698, -2)) #psuedo
gam1[2,2] ~ dnorm(-1.89 , pow(1.4825, -2))#psuedo
gam1[3,2] ~ dnorm(2.52 , pow(0.4313, -2))#psuedo

# second model
gam2[1,1] ~ dnorm(2.12 , pow(0.4664, -2)) # psuedo
gam2[2,1] ~ dnorm(-0.878 , pow(2.5838, -2)) # psuedo
gam2[3,1] ~ dnorm(2.57 , pow(0.4632, -2)) #psuedo
gam2[1,2] ~ dnorm(2.12 , pow(0.4664, -2))
gam2[2,2] ~ dnorm(-0.878 , pow(2.5838, -2))
gam2[3,2] ~ dnorm(2.57 , pow(0.4632, -2))


beta <- equals(m, 1) * beta1[,1] + equals(m, 2) * beta1[,1]

# first model
beta1[1,1] ~ dnorm(0.73, pow(0.437, -2))
beta1[2,1] ~ dnorm(3.06, pow(1.8636, -2))
beta1[3,1] ~ dnorm(-0.826, pow(0.324, -2))
beta1[1,2] ~ dnorm(0.73, pow(0.437, -2))  # psuedo
beta1[2,2] ~ dnorm(3.06, pow(1.8636, -2)) # psuedo
beta1[3,2] ~ dnorm(-0.826, pow(0.324, -2))# psuedo


# second model
beta2[1,1] ~ dnorm(0.73, pow(0.437, -2)) # psuedo
beta2[2,1] ~ dnorm(3.06, pow(1.8636, -2))   # psuedo
beta2[3,1] ~ dnorm(-0.826, pow(0.324, -2))# psuedo
beta2[1,2] ~ dnorm(0.727, pow(0.4371, -2))
beta2[2,2] ~ dnorm(6.40, pow(2.415, -2))
beta2[3,2] ~ dnorm(-0.828, pow(0.3361, -2))

phi <- equals(m, 1) * phi1[,,1] + equals(m, 2) *phi1[,,1]


phi1[1,1,1] ~ dnorm(0.167, pow(0.450, -2))
phi1[2,1,1] ~ dnorm(3.40, pow(2.449, -2))
phi1[3,1,1] ~ dnorm(3, pow(0.541, -2))
phi1[1,2,1] ~ dnorm(0.889, pow(0.5219, -2))
phi1[2,2,1] ~ dnorm(1.55, pow(1.5969, -2))
phi1[3,2,1] ~ dnorm(-2.26, pow(0.5698, -2))
phi1[1,3,1] ~ dnorm(2.48, pow(0.4249, -2))
phi1[2,3,1] ~ dnorm(0.116, pow(1.6815, -2))
phi1[3,3,1] ~ dnorm(0.0592, pow(0.4539, -2))
phi1[1,1,2] ~ dnorm(0.167, pow(0.450, -2)) # psuedo
phi1[2,1,2] ~ dnorm(3.40, pow(2.449, -2)) # psuedo
phi1[3,1,2] ~ dnorm(3, pow(0.541, -2)) # psuedo
phi1[1,2,2] ~ dnorm(0.889, pow(0.5219, -2)) # psuedo
phi1[2,2,2] ~ dnorm(1.55, pow(1.5969, -2)) # psuedo
phi1[3,2,2] ~ dnorm(-2.26, pow(0.5698, -2)) # psuedo
phi1[1,3,2] ~ dnorm(2.48, pow(0.4249, -2)) # psuedo
phi1[2,3,2] ~ dnorm(0.116, pow(1.6815, -2)) # psuedo
phi1[3,3,2] ~ dnorm(0.0592, pow(0.4539, -2)) # psuedo




# second model
phi2[1,1,1] ~ dnorm(2.04, pow(0.1725, -2)) # psuedo
phi2[2,1,1] ~ dnorm(3.40, pow(2.449, -2)) #psuedo
phi2[3,1,1] ~ dnorm(3, pow(0.541, -2)) #psuedo
phi2[1,2,1] ~ dnorm(0.889, pow(0.5219, -2)) #psuedo
phi2[2,2,1] ~ dnorm(3.59, pow(2.1774, -2)) # psuedo
phi2[3,2,1] ~ dnorm(-2.26, pow(0.5698, -2)) # psuedo
phi2[1,3,1] ~ dnorm(2.48, pow(0.4249, -2)) # psuedo
phi2[2,3,1] ~ dnorm(0.116, pow(1.6815, -2)) # psuedo
phi2[3,3,1] ~ dnorm(0.85, pow(0.1290, -2)) # psuedo
phi2[1,1,2] ~ dnorm(2.04, pow(0.1725, -2))
phi2[2,1,2] ~ dunif(-4, 4)
phi2[3,1,2] ~ dunif(-4, 4)
phi2[1,2,2] ~ dunif(-4, 4)
phi2[2,2,2] ~ dnorm(3.59, pow(2.1774, -2))
phi2[3,2,2] ~ dunif(-4, 4)
phi2[1,3,2] ~ dunif(-4, 4)
phi2[2,3,2] ~ dunif(-4, 4)
phi2[3,3,2] ~ dnorm(0.85, pow(0.1290, -2))

  for(i in 1:(nspec)){ 
    p[i] <- (exp(lp[i])) / (1 + exp(lp[i])) #This is the spot to fix correct?
  }

lp <- equals(m, 1) * lp1[,1] + equals(m, 2) * lp1[,1]

# first model
lp1[1,1] ~ dnorm(-1.009, pow(0.022, -2))
lp1[2,1] ~ dnorm(-3.6154, pow( 0.081, -2))
lp1[3,1] ~ dnorm(-1.337, pow(0.026, -2))
lp1[1,2] ~ dnorm(-1.009, pow(0.022, -2))# psuedo
lp1[2,2] ~ dnorm(-3.6154, pow( 0.081, -2))# psuedo
lp1[3,2] ~ dnorm(-1.337, pow(0.026, -2))# psuedo



# second model

lp2[1,1] ~ dnorm(-1.009, pow(0.0218, -2)) # psuedo
lp2[2,1] ~ dnorm(-3.6154, pow(0.0828, -2)) # psuedo
lp2[3,1] ~ dnorm(-1.337, pow(0.026, -2)) # psuedo
lp2[1,2] ~ dnorm(-1.009, pow(0.0218, -2))
lp2[2,2] ~ dnorm(-3.6154, pow(0.0828, -2))
lp2[3,2] ~ dnorm(-1.337, pow(0.026, -2))

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

