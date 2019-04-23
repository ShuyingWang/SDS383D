library(mosaic)
library(lme4)
library(tidyverse)
library(lattice)
library(MASS)
library(MCMCpack)

cheese = read.csv("/Users/zhendongwang/Desktop/R/cheese.csv")

## Linear regression
hlm3 = lmer(log(vol) ~ (1 + log(price) + disp + log(price):disp | store), data=cheese)
# Prior hyper-parameter for mu
nu = colMeans(coef(hlm3)$store)[c(4, 1, 2, 3)]
Kappa = 100 * diag(4)

## Build X and Y matrix
cheese$store = factor(cheese$store)
levels(cheese$store) <- c(1:88)
cheese = cheese[order(cheese$store), ]
M = as.numeric(summary(cheese$store))
n = 88
N = sum(M)
cM = cumsum(M)
mM = max(M)
p = log(cheese$price)
y = log(cheese$vol)
d = cheese$disp
l = cheese$store
k = 4

Y = matrix(nrow = n, ncol = mM)
X = array(dim = c(n, mM, k))
for (i in 1:n){
  Y[i, ] = c(y[(cM[i]-M[i]+1):cM[i]], rep(0, mM-M[i]))
  X[i, , 1] = c(rep(1, M[i]), rep(0, mM-M[i])) 
  X[i, , 2] = c(p[(cM[i]-M[i]+1):cM[i]], rep(0, mM-M[i]))
  X[i, , 3] = c(d[(cM[i]-M[i]+1):cM[i]], rep(0, mM-M[i]))
  X[i, , 4] = c(d[(cM[i]-M[i]+1):cM[i]] * p[(cM[i]-M[i]+1):cM[i]], rep(0, mM-M[i]))
}
Sxx = array(rep(0, n*k*k), dim = c(n, k, k))
Sxy = matrix(rep(0, n*k), nrow = n, ncol = k)
for (i in 1:n){
  for (j in 1:M[i]){
    Sxx[i, , ] = Sxx[i, , ] + X[i, j, ] %*% t(X[i, j, ])
    Sxy[i, ] = Sxy[i, ] + Y[i, j] * X[i, j, ]
  }
}

## Set initial values
mu = nu
lambda = 1
Sigma = matrix(c(4.99925, -0.94*2.2359*2.1732, 0.48*2.2359*0.9807, -0.33*2.2359*0.8336, 
                 -0.94*2.2359*2.1732, 4.72274, -0.56*2.1732*0.9807, 0.39*2.1732*0.8336, 
                 0.48*2.2359*0.9807, -0.56*2.1732*0.9807, 0.96184, -0.97*0.9807*0.8336, 
                 -0.33*2.2359*0.8336, 0.39*2.1732*0.8336, -0.97*0.9807*0.8336, 0.69495), nrow = 4, ncol = 4)
Beta = matrix(nrow = n, ncol = k)

## Gibbs sampling
burn = 5000
for (t in 1:burn){
  print(t)
  IS = solve(Sigma)
  b = 0
  
  # Update Beta
  for(i in 1:n){
    V = solve(IS + lambda * Sxx[i, , ])
    Mu = V %*% (lambda * Sxy[i, ] + IS %*% mu)
    Beta[i, ] = mvrnorm(1, Mu, V)
    
    b = b + sum((Y[i, 1:M[i]] - X[i, 1:M[i], ]%*%Beta[i, ])^2)
  }
  # Update lambda
  lambda = rgamma(1, (N+1)/2, b/2 + 1/2)
  
  # Update mu
  V = solve(n * IS + solve(Kappa))
  Mu = V %*% (solve(Kappa) %*% nu + IS %*% apply(Beta, 2, sum))
  mu = mvrnorm(1, Mu, V)
  
  # Update Sigma
  Mu = matrix(rep(mu,n), byrow=TRUE, nrow=n, ncol=k)
  Q = t(Beta - Mu) %*% (Beta - Mu)
  Sigma = riwish(n + 1, diag(k) + Q)
}

Sigma_sample = array(dim = c(B, k, k))
sigma_sample = c()
mu_sample = matrix(nrow = B, ncol = k)
Beta_sample = array(dim = c(B, n, k))

B = 10000
for (t in 1:B){
  print(t)
  IS = solve(Sigma)
  b = 0
  
  # Update Beta
  for(i in 1:n){
    V = solve(IS + lambda * Sxx[i, , ])
    Mu = V %*% (lambda * Sxy[i, ] + IS %*% mu)
    Beta[i, ] = mvrnorm(1, Mu, V)
    
    b = b + sum((Y[i, 1:M[i]] - X[i, 1:M[i], ]%*%Beta[i, ])^2)
  }
  # Update lambda
  lambda = rgamma(1, (N+1)/2, b/2 + 1/2)
  
  # Update mu
  V = solve(n * IS + solve(Kappa))
  Mu = V %*% (solve(Kappa) %*% nu + IS %*% apply(Beta, 2, sum))
  mu = mvrnorm(1, Mu, V)
  
  # Update Sigma
  Mu = matrix(rep(mu,n), byrow=TRUE, nrow=n, ncol=k)
  Q = t(Beta - Mu) %*% (Beta - Mu)
  Sigma = riwish(n + 1, diag(k) + Q)
  
  # Record samples
  Sigma_sample[t, , ] = Sigma
  sigma_sample[t] = lambda^(-0.5)
  mu_sample[t, ] = mu
  Beta_sample[t, , ] = Beta
}

## Results
sigma_hat = mean(sigma_sample)
Beta_hat = matrix(nrow = n, ncol = k)
error = 0
for (i in 1:n){
  Beta_hat[i, ] = colMeans(Beta_sample[ , i, ])
  error = error + sum((Y[i, 1:M[i]] - X[i, 1:M[i], ]%*%Beta[i, ])^2)
}
error = error / N


