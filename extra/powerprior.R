#############################
# NCC Sim Paper
# Power prior - Toy examples
# December 2024
#############################

rm(list=ls())

# install.packages("NPP")
library(NPP) 

# Consider three scenarios, where we vary the mean differences between concurrent and non-concurrent/historical data

# Scenario 1: current and historical means are similar
mean1_scenario1 <- 10 # mean of current data
mean0_scenario1 <- 10 # mean of historical data

# Scenario 2: current mean is larger than historical mean
mean1_scenario2 <- 10 # mean of current data
mean0_scenario2 <- 8 # mean of historical data

# Scenario 3: current mean is smaller than historical mean
mean1_scenario3 <- 10 # mean of current data
mean0_scenario3 <- 12 # mean of historical data

# Sample sizes
n1 <- 100 # Sample size for current data
n0 <- n1*2 # Sample size for historical data
var1 <- 1 # Variance of current data
var0 <- 1 # Variance of historical data

# Function to run the model with different delta priors and mean scenarios
run_model <- function(mean1, mean0, delta_alpha, delta_beta) {
 
 # Simulate data 
 data.Cur <- rnorm(n1, mean = mean1, sd = sqrt(var1)) # current data
 data.Hist <- rnorm(n0, mean = mean0, sd = sqrt(var0)) # historical data
 
 # Define the statistics
 CompStat <- list(
 n1 = n1, mean1 = mean1, var1 = var1,
 n0 = n0, mean0 = mean0, var0 = var0
 )
 
 # Power prior settings:
 # prior for the mean of the control response: we assume a normal prior with parameters mu0 (mean) and tau2 (variance)
 # prior for variance of the control response: we assume an inverse-gamma prior, with parameters alpha and beta
 # prior for delta: this refers to the weight of the power prior approach. We consider a beta prior with parameters delta.alpha and delta.beta
 
 # Define prior settings
 prior <- list(
 a = 1.5, # Power prior for variance
 mu0 = 0, # Prior mean for mu
 tau2 = 10, # Prior variance for mu
 alpha = 2, # Shape parameter for variance prior
 beta = 2, # Scale parameter for variance prior
 delta.alpha = delta_alpha, 
 delta.beta = delta_beta 
 )
 
 # Run the MCMC simulation  
 results <- NormalNPP_MCMC(
 Data.Cur = data.Cur,
 Data.Hist = data.Hist,
 CompStat = CompStat,
 prior = prior,
 MCMCmethod = 'RW', # Default: Random Walk method
 rw.logit.delta = 1, # Default: variance for random walk if delta is sampled
 nsample = 10000, # Default: Total MCMC iterations
 control.mcmc = list(
 delta.ini = NULL, # Default: Initial delta
 burnin = 2000, # Default: Burn-in iterations
 thin = 2 # Default: Thinning interval
 )
 )
 
 return(results)
}

# Setting 1: Very tight prior for delta -> Prior for full borrowing
#############
results_scenario1 <- run_model(mean1_scenario1, mean0_scenario1, 9000, 1000)
results_scenario2 <- run_model(mean1_scenario2, mean0_scenario2, 9000, 1000)
results_scenario3 <- run_model(mean1_scenario3, mean0_scenario3, 9000, 1000)

# tight but around 0.5
# results_scenario1 <- run_model(mean1_scenario1, mean0_scenario1, 10000, 10000)
# results_scenario2 <- run_model(mean1_scenario2, mean0_scenario2, 10000, 10000)
# results_scenario3 <- run_model(mean1_scenario3, mean0_scenario3, 10000, 10000)

# Results setting 1
# Scenario 1 (Similar means for current and historical data)
summary(results_scenario1$mu) # Posterior summary of mu
summary(results_scenario1$delta) # Posterior summary of delta

# Scenario 2 (current mean larger than historical mean):
summary(results_scenario2$mu)
summary(results_scenario2$delta)

# Scenario 3 (current mean smaller than historical mean):
summary(results_scenario3$mu)
summary(results_scenario3$delta)

# Setting 2: Less tight prior for delta
#############
results_scenario1 <- run_model(mean1_scenario1, mean0_scenario1, 100, 100)
results_scenario2 <- run_model(mean1_scenario2, mean0_scenario2, 100, 100)
results_scenario3 <- run_model(mean1_scenario3, mean0_scenario3, 100, 100)

# Results setting 2
# Scenario 1 (Similar means for current and historical data):
summary(results_scenario1$mu) # Posterior summary of mu
summary(results_scenario1$delta) # Posterior summary of delta

# Scenario 2 (current mean larger than historical mean):
summary(results_scenario2$mu)
summary(results_scenario2$delta)

# Scenario 3 (current mean smaller than historical mean):
summary(results_scenario3$mu)
summary(results_scenario3$delta)

# Setting 3: Very wide prior for delta -> this is the default in the pkg
#############
results_scenario1 <- run_model(mean1_scenario1, mean0_scenario1, 1, 1)
results_scenario2 <- run_model(mean1_scenario2, mean0_scenario2, 1, 1)
results_scenario3 <- run_model(mean1_scenario3, mean0_scenario3, 1, 1)

# Results setting 3
# Scenario 1 (Similar means for current and historical data):
summary(results_scenario1$mu) # Posterior summary of mu
summary(results_scenario1$delta) # Posterior summary of delta

# Scenario 2 (current mean larger than historical mean):
summary(results_scenario2$mu)
summary(results_scenario2$delta)

# Scenario 3 (current mean smaller than historical mean):
summary(results_scenario3$mu)
summary(results_scenario3$delta)

# Setting 4: Prior for no borrowing
#############
results_scenario1 <- run_model(mean1_scenario1, mean0_scenario1, 0.1, 1000)
results_scenario2 <- run_model(mean1_scenario2, mean0_scenario2, 0.1, 1000)
results_scenario3 <- run_model(mean1_scenario3, mean0_scenario3, 0.1, 1000)

# Results setting 4
# Scenario 1 (Similar means for current and historical data):
summary(results_scenario1$mu) # Posterior summary of mu
summary(results_scenario1$delta) # Posterior summary of delta

# Scenario 2 (current mean larger than historical mean):
summary(results_scenario2$mu)
summary(results_scenario2$delta)

# Scenario 3 (current mean smaller than historical mean):
summary(results_scenario3$mu)
summary(results_scenario3$delta)
