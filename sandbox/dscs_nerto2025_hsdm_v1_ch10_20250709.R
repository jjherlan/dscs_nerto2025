library(unmarked)
library(jagsUI)
library(R2WinBUGS)
library(R2OpenBUGS)

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

# 10.4 A slightly more complex site-occupancy model with covariates

# Choose sample sizes and prepare obs. data array y
set.seed(1) # So we get all the same data
M <- 100 # number of sites
J <- 3 # number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Create a covariate called vegHt

vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience

# Choose parameter values for occupancy model and compute occupancy

beta0 <- 0 # logit-scale intercept
beta1 <- 3 # logit-scale slope for vegHt
psi <- plogis(beta0 + beta1 * vegHt) # occupancy probability
#plot(vegHt, psi, type = "l", ylim = c(0,1), lwd = 3, xlab = "Vegetation height", ylab = "Occupancy probability") 
# plot psi relationship

# Now visit each site and observe presence/absence perfectly

z <- rbinom(n = M, size = 1, prob = psi) # occupancy state (true presence/absence)

# Look at data so far

table(z)

# Plot the true system state

par(mfrow = c(1, 3), mar = c(5, 5, 2, 2), cex.axis = 1.5, cex.lab = 1.5)

plot(vegHt, z, xlab = "Vegetation Height", ylab = "True presence/absence (z)", frame = F, 
     main = "True Occupancy State", pch = 19, cex = 1.5, col = ifelse(z == 1, "blue", "red"))

plot(function (x) plogis(beta0 + beta1*x), -1, 1, add = T, lwd = 3, col = "red")

# Create a covariate called wind

wind <- array(runif(M*J, -1, 1), dim = c(M, J)) # random wind covariate

# Choose parameter values for detection model and compute detection probability

alpha0 <- -2 # logit-scale intercept
alpha1 <- -3 # logit-scale slope for wind

p <- plogis(alpha0 + alpha1*wind) # detection probability

# plot(p ~ wind, ylim = c(0,1)) # plot detection probability relationship

# Take J = 3 presence/ansence measurements at each site

for (j in 1:J) {
  y[,j] <- rbinom(M, z, p[,j]) # detection state (presence/absence)
}

sum(apply(y, 1, max)) # observed number of occupied sites with at least one detection)

# Plot observed data and true effect of wind on detect probability

plot(wind, y, xlab="Wind", ylab="Observed det./nondetection data (y)", frame = F, ces = 1.5)
plot(function (x) plogis(alpha0 + alpha1*x), -1, 1, add = T, lwd = 3, col = "red")

# Look at the data: occupancy, tru presence/absence (z), and measurements (y)

cbind(psi = round(psi, 2), z = z, y1 = y[,1], y2 = y[,2], y3 = y[,3]) # first 6 sites (rows) of the data set

# Create factors
time <- matrix(rep(as.character(1:J), M), ncol = J, byrow = T) # time factor
hab <- c(rep("A", 33), rep("B", 33), rep("C", 34)) # Must have M = 100

# Load unmarked, format data and summarize

library(unmarked)

umf <- unmarkedFrameOccu(y = y, # Pres/Abs meausrements
                         siteCovs = data.frame(vegHt = vegHt, hab = hab), # site-specific covariates
                         obsCovs = list(wind = wind, time = time)) # obs-specific covariates

summary(umf) # Summarize the data

# Fit the model extract estimates
# Detection covariates follow first tilde, then occupancy covariates

summary(fm1.occ <- occu(~ wind ~ vegHt, data = umf)) # Fit the model

# Predict occupancy as function of covs (with 95% CIs)

# Add truth from data simulation (below for full code to produce fig. 10-2)

newdat <- data.frame(vegHt = seq(-1, 1, 0.01))

pred.occ <- predict (fm1.occ, type = "state", newdata = newdat)

newdat <- data.frame(wind = seq(-1, 1, 0.1))

pred.det <- predict(fm1.occ, type = "det", newdata = newdat)

# Predictions for specified values of vegHt, say 0.2 and 2.1

newdat <- data.frame(vegHt = c(0.2, 2.1))

predict(fm1.occ, type = "state", newdata = newdat, append = T)

# for values of wind of -1 to 1

newdat <- data.frame(wind=seq(-1, 1, 0.5))

predict(fm1.occ, type ="det", newdata = newdat, append = T)

# Fit detection-naive GLM to observed occurrence and plot comparison

summary(fm.glm <- glm(apply(y, 1, max) ~ vegHt, family = binomial))

plot(vegHt, apply(y, 1, max), xlab = "Vegetation height", 
     ylab = "Observed occurrence ('ever observed?;)", frame = F, cex = 1.5)

plot(function(x) plogis(beta0 + beta1*x), -1, 1, add = T, lwd = 3, col = "red")

lines(vegHt, predict(fm.glm,,"response"), type = "l", lwd = 3)

lines(vegHt, predict(fm1.occ, type = "state")[,1], col="blue", lwd = 3)

legend(-1, 0.9, c("Truth", "'LR' with p", "LR without p"), col = c("red", "blue", "black"), lty = 1, lwd = 3, cex = 1.2)

ranef(fm1.occ)

(psi1 <- predict(fm1.occ, type="state")[1,1])

(p1 <- predict(fm1.occ, type="det")[c(1:3),1])

(z1 <- (psi1*prod(1 - p1)) / ((1 - psi1) + psi1*prod(1 - p1)))

# Define the function for finite-sample number and proportion of occupied sites

fs.fn <- function(fm){
  Nocc <- sum(ranef(fm)@post[,2,])
  psi.fs <- Nocc / nrow(fm@data@y)
  out <- c(Nocc = Nocc, psi.fs = psi.fs)
  return(out)
}

# Bootstrap the function
fs.hat <- fs.fn(fm1.occ)           # Point estimate
pb.fs <- parboot(fm1.occ, fs.fn, nsim=10000, report=2) # Takes a while (33 min)
system.time(pb.fs <- parboot(fm1.occ, fs.fn, nsim=100, report=10)) # quicker

# Predictions for specified values of vegHt, say 0.2 and 2.1
newdat <- data.frame(vegHt=c(0.2, 2.1))
predict(fm.occ1, type="state", newdata=newdat, append = TRUE)


# Summarize bootstrap distributions
summary(pb.fs@t.star)




































































 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 