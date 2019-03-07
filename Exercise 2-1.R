## The conjugate Gaussian linear model (D)

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
m = c(0.01177, 0.20651)
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

