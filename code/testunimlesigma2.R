## Test UniMLESigma.R
rm(list = ls())
set.seed(1)
source('UniMLESigma.R')
n <- 1000
p <- 0.5
S <- rbinom(n, 1, p)

# y = b0 + b1*x + (1 + g1*x) N(0, 2)

b0 <- 1
b1 <- 1
g1 <- 0.5
sigma <- 2
x <- runif(n, 0, 1)
y <- rep(0, n)
for (i in 1:n){
  if (S[i] == 1) {
    y[i] <- b0 + b1*x[i]+ (1+g1*x[i])* rnorm(1, 0, sigma)
  } else {
    y[i] <- b0 + b1*x[i]+ (1+g1*x[i])* rnorm(1, 0, sigma)
  }
}

## mean regression of y on x
print(summary(lm(y ~ x)))
plot(y ~ x)

mod5 <- QRGradient(y, S, x, tau = 0.5)
abline(mod5$param[1:2])
print(mod5$param)

mod7 <- QRGradient(y, S, x, tau = 0.7)
abline(mod7$param[1:2])
print(mod7$param)

mod9 <- QRGradient(y, S, x, tau = 0.9)
abline(mod9$param[1:2])
print(mod9$param)

mod3 <- QRGradient(y, S, x, tau = 0.3)
abline(mod3$param[1:2])
print(mod3$param)

mod1 <- QRGradient(y, S, x, tau = 0.1)
abline(mod1$param[1:2])
print(mod1$param)
