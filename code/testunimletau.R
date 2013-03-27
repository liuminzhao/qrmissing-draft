## Test UniMLETau.R
rm(list = ls())
set.seed(1)
source('UniMLETau.R')
n <- 200
p <- 0.5
S <- rbinom(n, 1, p)
b01 <- 1
b00 <- -1
b11 <- 1
b10 <- -1
sigma1 <- 1
sigma0 <- 1
x <- runif(n, 0, 2)
y <- rep(0, n)
for (i in 1:n){
  if (S[i] == 1) {
    y[i] <- rnorm(1, b01 + x[i] * b11, sigma1)
  } else {
    y[i] <- rnorm(1, b00 + x[i] * b10, sigma0)
  }
}
y1 <- y[S == 1]
y0 <- y[S == 0]

## mean regression of y on x
summary(lm(y ~ x))
plot(y ~ x)

source('UniMLETau.R')
mod <- QRGradient(y, S, x, tau = 0.5, niter = 1000)

mod7 <- QRGradient(y, S, x, tau = 0.7)

mod9 <- QRGradient(y, S, x, tau = 0.9)

mod3 <- QRGradient(y, S, x, tau = 0.3)

mod2 <- QRGradient(y, S, x, tau = 0.2)
mod1 <- QRGradient(y, S, x, tau = 0.1, niter = 1000)

Rprof("Rprof.out")

abline(mod1$gamma)
abline(mod3$gamma)
abline(mod7$gamma)
abline(mod9$gamma)

Rprof(NULL)
summaryRprof()

myplot(mod)

