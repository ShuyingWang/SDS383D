## Polls model 1 ##

library(mvtnorm)
library(truncnorm)
library(MCMCpack)
library(fastDummies)

## Load in data
polls = read.csv("/Users/zhendongwang/Desktop/R/polls.csv")
polls <- polls[complete.cases(polls), ]

polls$state = factor(polls$state)
levels(polls$state) <- c(1:49)
polls = polls[order(polls$state), ]   # Sort data by state

M = as.numeric(summary(polls$state))   # sample size for each state
cM = cumsum(M)
n = length(M)   # number of states
N = sum(M)   # Whole sample size
k = 11 # dim of X

polls = dummy_cols(polls, select_columns = c('edu', 'age'))   # Code dummy variables
polls$weight = (polls$weight - mean(polls$weight)) / sd(polls$weight)   # Normalize weight

Y = polls$bush
X = cbind(polls[, 8], polls[, 9], polls[, 10], polls[, 11], polls[, 12], polls[, 13], polls[, 14], polls[, 15], polls[, 16], polls[, 17], polls[, 18])
Sxx = array(rep(0, n*k*k), dim = c(n, k, k))
for (i in 1:n){
  for (j in M[i]:1){
    Sxx[i, , ] = Sxx[i, , ] + X[(cM[i]-j+1), ] %*% t(X[(cM[i]-j+1), ])
  }
}


## Function for gibbs samplers
Sampler_dp <- function(mu, Beta){
  dp = matrix(nrow = N, ncol = 2)
  for (i in 1:N){
    state = polls$state[i]
    z = mu[state] + sum(X[i, ] * Beta[state, ])
    dp[i, 2] = pnorm(z, 0, 1)
    dp[i, 1] = ifelse(Y[i] == 1,  rtruncnorm(1, a = 0, b = Inf, mean = z, sd = 1), 
                  rtruncnorm(1, a = -Inf, b = 0, mean = z, sd = 1))
  }
  return(dp)
}

Sampler_mu <- function(d, Beta, nu, lambda){
  mu = c()
  for (i in 1:n){
    mu[i] = rnorm(1, (lambda*nu + sum(d[(cM[i]-M[i]+1):cM[i]] - X[(cM[i]-M[i]+1):cM[i], ]%*%Beta[i, ])) / (M[i]+lambda), (M[i]+lambda)^(-0.5))
  }
  return(mu)
}

Sampler_Beta <- function(d, mu, phi, Sigma){
  Beta = matrix(nrow = n, ncol = k)
  IS = solve(Sigma)
  for (i in 1:n){
    V = solve(IS + Sxx[i, , ])
    A = IS %*% t(phi) + apply((d[(cM[i]-M[i]+1):cM[i]] - mu[i]) * X[(cM[i]-M[i]+1):cM[i], ], 2 , sum)
    Beta[i, ] = rmvnorm(1, V %*% A, V)
  }
  return(Beta)
}

Sampler_nu <- function(lambda, mu){
  nu = rnorm(1, mean(mu), (n*lambda)^(-0.5))
  return(nu)
}

Sampler_lambda <- function(mu, nu){
  lambda = rgamma(1, (n+1)/2, 1/2 + 1/2 * sum((mu - nu)^2))
  return(lambda)
}

Sampler_phi <- function(Sigma, Beta){
  phi = rmvnorm(1, colMeans(Beta), Sigma / n)
  return(phi)
}

Sampler_Sigma <- function(Beta, phi){
  Phi = matrix(rep(phi,n), byrow=TRUE, nrow=n, ncol=k)
  Q = t(Beta - Phi) %*% (Beta - Phi)
  Sigma = riwish(n + 1, diag(k) + Q)
  return(Sigma)
}

## initial values
mu = rnorm(n, 0, 1)
Beta = matrix(rnorm(k*n, 0, 1), nrow = n, ncol = k)
nu = 0
lambda = 1
Sigma = diag(k)
phi = rmvnorm(1, colMeans(Beta), Sigma / n)


## Gibbs sampling
burn = 5000
for (t in 1:burn){
  print(t)
  d  = Sampler_dp(mu, Beta)[ , 1]
  mu = Sampler_mu(d, Beta, nu, lambda)
  Beta = Sampler_Beta(d, mu, phi, Sigma)
  nu = Sampler_nu(lambda, mu)
  lambda = Sampler_lambda(mu, nu)
  phi = Sampler_phi(Sigma, Beta)
  Sigma = Sampler_Sigma(Beta, phi)
}

B = 10000
p_sample = matrix(nrow = B, ncol = N)
d_sample = matrix(nrow = B, ncol = N)
mu_sample = matrix(nrow = B, ncol = n)
Beta_sample = array(dim = c(B, n, k))

for (t in 1:B){
  print(t)
  dp  = Sampler_dp(mu, Beta)
  d = dp[ , 1]
  mu = Sampler_mu(d, Beta, nu, lambda)
  Beta = Sampler_Beta(d, mu, phi, Sigma)
  nu = Sampler_nu(lambda, mu)
  lambda = Sampler_lambda(mu, nu)
  phi = Sampler_phi(Sigma, Beta)
  Sigma = Sampler_Sigma(Beta, phi)
  
  p_sample[t, ] = dp[ , 2]
  d_sample[t, ] = d
  mu_sample[t, ] = mu
  Beta_sample[t, , ] = Beta
}

## Results
p_hat = colMeans(p_sample)
y_hat = ifelse(p_hat > 0.5, 1, 0)
error = sum(abs(Y - y_hat)) / N
hist(p_hat)





## Polls model 2 ##

library(mvtnorm)
library(truncnorm)
library(fastDummies)

## Load in data
polls = read.csv("/Users/zhendongwang/Desktop/R/polls.csv")
polls <- polls[complete.cases(polls), ]

polls$state = factor(polls$state)
levels(polls$state) <- c(1:49)
polls = polls[order(polls$state), ]   # Sort data by state

M = as.numeric(summary(polls$state))   # sample size for each state
cM = cumsum(M)
n = length(M)   # number of states
N = sum(M)   # Whole sample size
k = 10 # dim of X

polls = dummy_cols(polls, select_columns = c('edu', 'age'), remove_first_dummy = TRUE)   # Code dummy variables
polls$weight = (polls$weight - mean(polls$weight)) / sd(polls$weight)   # Normalize weight

Y = polls$bush
X = cbind(rep(1, N), polls[, 8], polls[, 9], polls[, 10], polls[, 11], polls[, 12], polls[, 13], polls[, 14], polls[, 15], polls[, 16])
Sxx = array(rep(0, n*k*k), dim = c(n, k, k))
for (i in 1:n){
  for (j in M[i]:1){
    Sxx[i, , ] = Sxx[i, , ] + X[(cM[i]-j+1), ] %*% t(X[(cM[i]-j+1), ])
  }
}

## Function for gibbs samplers
Sampler_zp <- function(Beta){
  zp = matrix(nrow = N, ncol = 2)
  for (i in 1:N){
    state = polls$state[i]
    xb = sum(X[i, ] * Beta[state, ])
    zp[i, 2] = pnorm(xb, 0, 1)
    zp[i, 1] = ifelse(Y[i] == 1,  rtruncnorm(1, a = 0, b = Inf, mean = xb, sd = 1), 
                      rtruncnorm(1, a = -Inf, b = 0, mean = xb, sd = 1))
  }
  return(zp)
}

Sampler_Beta <- function(z, phi, lambda){
  Beta = matrix(nrow = n, ncol = k)
  for (i in 1:n){
    V = solve(lambda * diag(k) + Sxx[i, , ])
    A = lambda*phi + apply(z[(cM[i]-M[i]+1):cM[i]] * X[(cM[i]-M[i]+1):cM[i], ], 2 , sum)
    Beta[i, ] = rmvnorm(1, V %*% t(A), V)
  }
  return(Beta)
}

Sampler_phi <- function(lambda, Beta){
  phi = rmvnorm(1, colMeans(Beta), diag(k) / (n*lambda))
  return(phi)
}

Sampler_lambda <- function(Beta, phi){
  Phi = matrix(rep(phi,n), byrow=TRUE, nrow=n, ncol=k)
  lambda = rgamma(1, (n*k-1)/2, 1/2 + 1/2 * sum((Beta - Phi)^2))
  return(lambda)
}


## initial values
Beta = matrix(rnorm(k*n, 0, 1), nrow = n, ncol = k)
lambda = 1
phi = rmvnorm(1, colMeans(Beta), diag(k) / (n*lambda))


## Gibbs sampling
burn = 5000
for (t in 1:burn){
  print(t)
  z  = Sampler_zp(Beta)[ , 1]
  Beta = Sampler_Beta(z, phi, lambda)
  phi = Sampler_phi(lambda, Beta)
  lambda = Sampler_lambda(Beta, phi)
}

B = 10000
p_sample = matrix(nrow = B, ncol = N)
Beta_sample = array(dim = c(B, n, k))
phi_sample = matrix(nrow = B, ncol = k)

for (t in 1:B){
  print(t)
  zp  = Sampler_zp(Beta)
  z = zp[ , 1]
  Beta = Sampler_Beta(z, phi, lambda)
  phi = Sampler_phi(lambda, Beta)
  lambda = Sampler_lambda(Beta, phi)
  
  p_sample[t, ] = zp[ , 2]
  phi_sample[t, ] = phi
  Beta_sample[t, , ] = Beta
}

## Results
p_hat = colMeans(p_sample)
y_hat = ifelse(p_hat > 0.5, 1, 0)
error = sum(abs(Y - y_hat)) / N
hist(p_hat)

