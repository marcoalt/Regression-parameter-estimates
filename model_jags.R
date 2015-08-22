model {
  for (i in 1:N)
  {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b1 * x[i] + b0
  }
  
  b1 ~ dnorm(0, 0.0001)
  b0 ~ dnorm(0, 0.0001)
  
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 10000)
}