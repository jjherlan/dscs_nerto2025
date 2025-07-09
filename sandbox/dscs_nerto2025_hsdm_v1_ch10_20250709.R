library(unmarked)
library(jagsUI)
library(R2WinBUGS)

# library(mgcv) # Not used in this code snippet, but can be used for GAMs

# 10.3 Simulation and analysis of the simlest possible site-occupancy model

# Load required libraries

# choose sample sizes and prepare observed array y

set.seed(24) # So we get all the same data
M <- 100 # number of sites
J <- 2 # number of presenc/absence measurements
y <- matrix(0, nrow = M, ncol = J) # to contain the obs. data
 
# Parameter values
 
psi <-0.8 # probability of occupancy or presence
p <- 0.5 # Probability of detection given occupancy
 
# Generate the presence/absence data (the truth)
z <- rbinom(n = M, size = 1, prob = psi) # R has no Bernoulli [occupancy state]
 
# Generate the detection nondetection data (i.e., presence/absence measurements)
 for (j in 1:J) {
   y[,j] <- rbinom(n = M, size = 1, prob = z*p) # R has no Bernoulli [detection state]
 }
 
# Look at data
sum(z) # true number of occupied sites
 
sum(apply(y, 1, max)) # observed number of occupied sites with at least one detection

# In this simulation, the species is present at 80% of the sites, 
# and we detect it at 50% of the sites where it is present. We detect it 
# at 61. The overall measurement error for the apparent number of occupied sites is thus 
# (86 - 61)/86 = 0.292, or -29.2%.Under the binomial model, we'd expect a combined detection 
# probability (over J surveys) of 1 - (1 - p)^J = 0.75, so the measurement error of -25%. 
# This difference between -29% and -25% is of course due to the sampling error inherent in the stochastic detection process. 

# Now we inspect our data set:

head(cbind(z = z, y1 = y[.1], y2 = y[,2])) # Truth and measurements for first 6 sites (rows) of the data set
 
# Now we can analyze the data using the unmarked package using the function 'occu',
# where the linear model for detection is specified before tha tfor occupancy

umf <- unmarkedFrameOccu(y = y) # Create unmarked data frame
summary(umf) # Summarize the data
(fm1 <- occu(~1~1, data = umf))

# Bundle data and summarize bundle data

str(win.data <- list(y = y, M = nrow(y), J = ncol(y)) )

# Specify model in BUGS language

sink("model.txt")
cat("
model {
  # Priors
  psi ~ dunif(0, 1) # Prior for occupancy probability
  p ~ dunif(0, 1) # Prior for detection probability
  
  # Likelihood
  for (i in 1:M) { # Loop over sites
    z[i] ~ dbern(psi) # State model
    for (j in 1:J) { # Loop over replicate surveys
      y[i,j] ~ dbern(z[i] * p) # Observation model (only JAGS!)
#      y[i,j] ~ dbern(mu[i]) # for WinBUGS define 'strawman'
    }
#   mu[i] <- z[i] * p # Only WinBUGS 
  }
}
",fill = TRUE)
sink()

# Initial values

zst <- apply(y, 1, max) # Avoid data/model/inits conflict
inits <- function() {
  list(z = zst)
}
 
# Parameters restored
params <- c("psi", "p")

# MCMC settings

ni <- 5000 ; nt <- 1 ; nb <- 1000 ; nc <- 3 # Number of iterations

# Call JAGS and summarize posteriors

fm2 <- jags(data = win.data, inits = inits, parameters.to.save = params,
          model.file = "model.txt", n.chains = nc, n.thin = nt, n.burnin = nb,
          n.iter = ni, DIC = TRUE)
print(fm2, digits = 3)







 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 