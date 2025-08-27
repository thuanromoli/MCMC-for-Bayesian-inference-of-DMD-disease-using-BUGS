library(nimble)
library(igraph)
library(coda)
library(R6)
load("DMDdata")

##########################
####### Model setup  #####
##########################

# Model specification
DMDCode1 <- nimbleCode({
  
  # Specify the likelihood:
  for (i in 1:N){
    Y[i] ~ dbin(pi[i],ni[i])
    #pi[i] <- expit(beta_0)
    #pi[i] <- expit(beta_0 + beta_hp * ( ( X_hp[i]-mean(X_hp[]) ) /sd(X_hp[]) ))
    #pi[i] <- expit(beta_0 + beta_hi * ( ( X_hi[i]-mean(X_hi[]) ) /sd(X_hi[]) ))
    pi[i] <- expit(beta_0 +
                     beta_hp * ( ( X_hp[i]-mean(X_hp[]) ) /sd(X_hp[]) ) +
                     beta_hi * ( ( X_hi[i]-mean(X_hi[]) ) /sd(X_hi[]) ))
  }
  # Prior specification: (N(0,tau) where tau = 1/sigma^2)
  beta_0 ~ dnorm(0,1/25)
  beta_hp ~ dnorm(0,1/25)
  beta_hi ~ dnorm(0,1/25)
  
})

# Constant values in the model
DMDConsts <- list(N=15)

# Data values
DMDData <- list(ni = DMDdata$ni,
                X_hp = DMDdata$X_hp,
                X_hi = DMDdata$X_hi,
                Y = DMDdata$ncarriers)

# Initial values before building the model                 
DMDInits <- list(beta_0 = 1, beta_hp = 1, beta_hi = 1)

# To build the model
DMDmod <- nimbleModel(code = DMDCode1, name = "DMDmod", constants = DMDConsts,
                      data = DMDData, inits = DMDInits)
# To compile the model
CDMDmod <- compileNimble(DMDmod)


# Set up the monitored quantities
DMDmodConfig <- configureMCMC(DMDmod,enableWAIC = TRUE, monitors = c('beta_0','beta_hp','beta_hi'), print = TRUE) 
# Build the MCMC algorithm
DMDmodMCMC <- buildMCMC(DMDmodConfig)
# Compile the MCMC chain 
CDMDmodMCMC <- compileNimble(DMDmodMCMC, project = DMDmod)

set.seed(10) # for replicability

DMDInits <- list(list(beta_0 = 0, beta_hp = 0, beta_hi = 0),
                 list(beta_0 = 7, beta_hp = 7, beta_hi = 7),
                 list(beta_0 = 1, beta_hp = 1, beta_hi = 1),
                 list(beta_0 = -10, beta_hp = -10, beta_hi = -10))

posterior <- runMCMC(CDMDmodMCMC, niter = 10000, thin=1, nburnin=1, 
                     summary = TRUE, samples = TRUE, nchains=4, 
                     samplesAsCodaMCMC=TRUE, inits = DMDInits) 

combinedchains <- mcmc.list(posterior$samples$chain1,
                            posterior$samples$chain2,
                            posterior$samples$chain3,
                            posterior$samples$chain4)

posterior$summary$all.chains

############################
#### Model diagnostic  #####
############################

# Checks for convergence and burn-in
plot(combinedchains)
gelman.plot(combinedchains)

# Checks for enough sample sizes
# ACF plots
a <- acf(posterior$samples$chain1[, "beta_0"]) # plot autocorrelation of alpha sample
b <- acf(posterior$samples$chain1[, "beta_hp"]) # plot autocorrelation of beta sample
c <- acf(posterior$samples$chain1[, "beta_hi"]) # plot autocorrelation of tau sample

par(mfrow = c(1, 3))
plot(a, main = "ACF of beta_0")
plot(b, main = "ACF of beta_hp")
plot(c, main = "ACF of beta_hi")

# Effective sample size
effectiveSize(combinedchains)
# Monte Carlo error
summary(combinedchains)

# checks for all betas if Time-Series SE / SD < 0.05
statistics <- summary(combinedchains)[[1]]
j <- c("beta_0","beta_hi","beta_hp")
for (i in 1:3) {
  print(paste(j[i],"has a ratio of",statistics[i,4] / statistics[i,2] ))
}



#calculateWAIC(CDMDmodMCMC)# offline approach to WAIC


############################
#### Poseterior checks  ####
############################

# loading data
ni = DMDdata$ni
X_hp = DMDdata$X_hp
X_hi = DMDdata$X_hi
Y = DMDdata$ncarriers

# Number of posterior predictive samples

nsamples <- 10000

# Initialising lists to store 10000 samples for the 15 different hospitals
Y_samples <-  vector(mode = "list", length = nsamples)

for (i in 1:nsamples) {
  # Generating 10000 posterior predictive samples
  # For each parameter, sample from the empirical distribution
  beta_0_sample <- sample(as.matrix(posterior$samples)[,1], 1)
  beta_hp_sample <- sample(as.matrix(posterior$samples)[,3], 1)
  beta_hi_sample <- sample(as.matrix(posterior$samples)[,2], 1)
  
  # Calculate probabilities for each hospital
  pi_sample <- expit(beta_0_sample +
                    beta_hp_sample * ( ( X_hp[]-mean(X_hp[]) ) /sd(X_hp[]) ) +
                    beta_hi_sample * ( ( X_hi[]-mean(X_hi[]) ) /sd(X_hi[]) ))
  
  # Generate a binomial sample for each hospital
  Y_samples[[i]] <- rbinom(15,ni,pi_sample)
}


# Initialising numeric vectors to store a distribution for the mean, median, and max carriers
Y_samples_mean <-  vector(mode = "numeric", length = nsamples)
Y_samples_median <-  vector(mode = "numeric", length = nsamples)
Y_samples_max <-  vector(mode = "numeric", length = nsamples)

# Calculates the mean, median, and max carriers for each sample
for (i in 1:nsamples) {
  Y_samples_mean[i] <- mean(Y_samples[[i]])
  Y_samples_median[i] <- median(Y_samples[[i]])
  Y_samples_max[i] <- max(Y_samples[[i]])
}

# Calculates the Bayesian p values
pval_mean <- min(length(which(Y_samples_mean>=mean(Y))) / nsamples,
                 length(which(Y_samples_mean<=mean(Y))) / nsamples)
pval_median <- min(length(which(Y_samples_median>=median(Y))) / nsamples,
                   length(which(Y_samples_median<=median(Y))) / nsamples)
pval_max <- min(length(which(Y_samples_max>=max(Y))) / nsamples,
                 length(which(Y_samples_max>=max(Y))) / nsamples)

print(paste("Bayesian p-value for the mean:", pval_mean))
print(paste("Bayesian p-value for the median:", pval_median))
print(paste("Bayesian p-value for the max:", pval_max))

par(mfrow = c(1, 3))
hist(Y_samples_mean, prob = TRUE, col = 'white', xlab = "Mean", main="Mean of carriers", breaks = 30)
abline(v=mean(Y),lty=2,col='red')
hist(Y_samples_median, prob = TRUE, col = 'white', xlab = "Median", main="Median of carriers",breaks = 30)
abline(v=median(Y),lty=2,col='red')
hist(Y_samples_max, prob = TRUE, col = 'white', xlab = "Max", main="Max of carriers",breaks = 30)
abline(v=max(Y),lty=2,col='red')

