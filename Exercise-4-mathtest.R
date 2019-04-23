library(mosaic)	# nice for its overloading of the mean() function

mathtest = read.csv('/Users/zhendongwang/Desktop/R/mathtest.csv')

# Ensure that school is treated as a categorical variable
mathtest$school = factor(mathtest$school)

# This won't work without mosaic loaded
schoolmeans = mean(mathscore~school, data=mathtest)
schoolsizes = as.numeric(summary(mathtest$school))

# Notice the extremes tend to be for schools with low sample sizes
plot(schoolsizes, schoolmeans, pch=19)

n = length(schoolsizes)
N = summary(mathtest$school)

y = list(c())
i = 1
for (k in 1:sum(N)){
  if (mathtest$school[k] == i) {
    y[[i]] = c(y[[i]], mathtest$mathscore[k])
  } else {
    i = i + 1
    y[[i]] = c(mathtest$mathscore[k])
  }
}

# Hyperparameters
m = mean(schoolmeans)
v = 10000
a = 0.1
b = 0.1
c = 0.1
d = 0.1

# Gibbs samplers
Sampler_mu <- function(theta, eta){
  omega = n*eta + 1/v
  mu_star = (eta * sum(theta) + m/v) / omega
  mu_new = rnorm(1, mu_star, omega^(-0.5))
  return(mu_new)
}
Sampler_lambda <- function(theta){
  a_star = a + sum(N)/2
  s = 0
  for (i in 1:n){
    s = s + sum((y[[i]] - theta[i])^2)
  }
  b_star = b + s/2
  lambda_new = rgamma(1, a_star, b_star)
  return(lambda_new)
}
Sampler_eta <- function(theta, mu){
  c_star = c + n/2
  d_star = d + 0.5 * sum((theta - mu)^2)
  eta_new = rgamma(1, c_star, d_star)
  return(eta_new)
}
Sampler_theta <- function(mu, lambda, eta){
  theta_new = c()
  for (i in 1:n){
    omega = N[i]*lambda + eta
    mu_star = (lambda * sum(y[[i]]) + eta * mu) / omega
    theta_new[i] = rnorm(1, mu_star, omega^(-0.5))
  }
  return(theta_new)
}

# Set initial values
mu = m
lambda = 0.01
eta = 0.01
theta = rep(m, n)

# 5000 iterations before burn in
for (t in 1:5000){
  mu = Sampler_mu(theta, eta)
  lambda = Sampler_lambda(theta)
  eta = Sampler_eta(theta, mu)
  theta = Sampler_theta(mu, lambda, eta)
}

# Vectors to contain samples
B = 10000
sample_mu = c()
sample_lambda = c()
sample_eta = c()
sample_theta = matrix(nrow = B, ncol = n)

# 10,000 iterations after burn in
for (t in 1:B){
  mu = Sampler_mu(theta, eta)
  lambda = Sampler_lambda(theta)
  eta = Sampler_eta(theta, mu)
  theta = Sampler_theta(mu, lambda, eta)
  
  sample_mu[t] = mu
  sample_lambda[t] = lambda
  sample_eta[t] = eta
  sample_theta[t, ] = theta
}

# Compute sample means
mu_hat = mean(sample_mu)
lambda_hat = mean(sample_lambda)
sigma_hat = lambda_hat^(-0.5)
eta_hat = mean(sample_eta)
tau_hat = (eta_hat * sigma_hat^2)^(-0.5)
theta_hat = colMeans(sample_theta)

error = c()
for (i in 1:n){
  error = c(error, (y[[i]] - theta_hat[i]) / sigma_hat)
}
hist(error, breaks = 30, freq = FALSE)
curve(expr = dnorm(x, 0, 1), add = TRUE)


# Shrinkage coefficient
k = abs(schoolmeans - theta_hat) / schoolmeans
plot(N, k, pch=19, cex = 0.7)
