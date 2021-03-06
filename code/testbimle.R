## Test BiMLE.R
## SIMULATE DATA

set.seed(1)
p <- 0.5
n <- 200
R <- rbinom(n, 1, p)

b01 <- 1
b00 <- -1
b11 <- 1
b10 <- -1
sigma1 <- 1
sigma0 <- 1

b02 <- 1
b12 <- -1
sigma2 <- 1

x <- runif(n, 0, 2)
y <- matrix(0, n, 2)
for (i in 1:n){
  if (R[i] == 1){
    y[i, 1] <- rnorm(1, b01 + x[i] * b11, sigma1)
    y[i, 2] <- rnorm(1, b02 + x[i] * b12, sigma2)
  } else {
    y[i, 1] <- rnorm(1, b00 + x[i] * b10, sigma0)
    y[i, 2] <- NA
  }
}

## MEAN REGRESSION AND PLOT
summary(lm(y[,1] ~ x))
plot(y[,1] ~ x)

summary(lm(y[,2] ~ x))
plot(y[,2] ~ x)

## options(warn = 2)
source('BiMLE.R')
mod <- BiQRGradient(y, R, x, tau = 0.4, niter = 1000)
print(mod$gamma1)
print(mod$gamma2)
print(mod$tau)
