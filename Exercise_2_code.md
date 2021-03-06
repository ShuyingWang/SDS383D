Exercise 2
================
Shuying Wang
3/7/2019

The conjugate Gaussian linear model (D)
---------------------------------------

``` r
# load in data
library(MASS)
gdpdata = read.csv("/Users/zhendongwang/Desktop/R/gdpgrowth.csv")
Y = gdpdata$GR6096
n = length(Y)
X = cbind(rep(1, n), gdpdata$DEF60)

# OLS
OLS <- lm(GR6096~DEF60, data = gdpdata)
lm1 <- function(x){
  y = 0.20651*x + 0.01177
  return(y)
}
plot(gdpdata[ ,"DEF60"], gdpdata[ ,"GR6096"], cex = 0.7)
curve(expr = lm1, add = TRUE)

# Bayesian linear model
# Prior hyper-parameters
m = c(0, 0)
K = diag(c(0.01, 0.01))
d = 1
eta = 1

# Posterior hyper-parameters
d_post = d + n
eta_post = eta + m%*%K%*%m + Y%*%Y - t(K%*%m + t(X)%*%Y) %*% solve(K + t(X)%*%X) %*% (K%*%m + t(X)%*%Y)
m_post = solve(K + t(X)%*%X) %*% (K%*%m +t(X)%*%Y)
K_post = K + t(X)%*%X

# Sample from posterior
N = 1000
omega_post = rgamma(N, d_post/2, rate = eta_post/2)
beta_post = matrix(nrow = N, ncol = 2)
for (i in 1:N){
  beta_post[i, ] = mvrnorm(1, m_post, solve(omega_post[i]*K_post))
}
beta = colMeans(beta_post)

# Graph
lm2 <- function(x){
  y = beta[2]*x + beta[1]
  return(y)
}
curve(expr = lm2, add = TRUE, col = 'red')
legend('topleft',lty=1,c("OLS",'Bayesian'),col=c('black','red'))
```

![](Exercise_2_code_files/figure-markdown_github/unnamed-chunk-1-1.png)

A heavy-tailed error model (C)
------------------------------

``` r
# load in data
library(MASS)
gdpdata = read.csv("/Users/zhendongwang/Desktop/R/gdpgrowth.csv")
Y = gdpdata$GR6096
n = length(Y)
X = cbind(rep(1, n), gdpdata$DEF60)

# Non-heavy-tailed model (result from The conjugate Gaussian linear model (D))
lm1 <- function(x){
  y = 0.1892664*x + 0.01197588
  return(y)
}
plot(gdpdata[ ,"DEF60"], gdpdata[ ,"GR6096"], cex = 0.7)
curve(expr = lm1, add = TRUE)

# Heavy-tailed error model
# Prior hyper-parameters
m = c(0, 0)
K = diag(c(0.001, 0.001))
d = 10
eta = 1
h = 1

# Functions for compute posterior hyper-parameters
Mean_beta <- function(delta){
  m = solve(t(X)%*%delta%*%X + K) %*% (t(X)%*%delta%*%Y + K%*%m)
  return(m)
}
Sigma_beta <- function(delta, omega){
  s = solve(omega * (t(X)%*%delta%*%X + K))
  return(s)
}
A_omega <- (d+n)/2 + 2
B_omega <- function(delta){
  b = (eta + Y%*%delta%*%Y + m%*%K%*%m - (t(t(X)%*%delta%*%Y + K%*%m) %*% solve(t(X)%*%delta%*%X + K) %*% (t(X)%*%delta%*%Y + K%*%m))) / 2
  return(b)
}
A_lambda <- (h+1)/2
B_lambda <- function(x, y, omega, beta){
  b = (h + omega * (y - x%*%beta)^2) / 2
  return(b)
}

# Gibbs Sampler
# Set initial values
beta = c(0, 0)
omega = 1
delta = diag(rep(0, n))

# First 1000 iterations before burn in
for (i in 1:1000){
  beta = mvrnorm(1, Mean_beta(delta), Sigma_beta(delta, omega))
  omega = rgamma(1, A_omega, B_omega(delta))
  for (j in 1:n){
    delta[j, j] = rgamma(1, A_lambda, B_lambda(X[j, ], Y[j], omega, beta))
  }
}

# 1000 iterations after burn in
Beta <- matrix(nrow = 1000, ncol = 2)
for (i in 1:1000){
  beta = mvrnorm(1, Mean_beta(delta), Sigma_beta(delta, omega))
  omega = rgamma(1, A_omega, B_omega(delta))
  for (j in 1:n){
    delta[j, j] = rgamma(1, A_lambda, B_lambda(X[j, ], Y[j], omega, beta))
  }
  Beta[i, ] = beta
}
beta_post = colMeans(Beta)

# Graph
lm2 <- function(x){
  y = beta_post[2]*x + beta_post[1]
  return(y)
}
curve(expr = lm2, add = TRUE, col = 'red')
legend('topleft',lty=1,c("Non-heavy-tailed",'heavy-tailed'),col=c('black','red'))
```

![](Exercise_2_code_files/figure-markdown_github/unnamed-chunk-2-1.png)
