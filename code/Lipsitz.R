Lipsitz <- function(X, y) {

}

alpha <- 0
n <- 30
p <- 0.5
R <- rbinom(n, 1, p)
x <- runif(n, 0, 2)
y <- matrix(0, n, 2)
for (i in 1:n){
    if (R[i] == 1){
        y[i, 1] <- 2 + x[i] +(1 + alpha*x[i])* rnorm(1)
        y[i, 2] <- 1 - x[i] - 0.5 * y[i, 1] + (1+alpha*x[i])* rnorm(1)
    } else {
        y[i, 1] <- -2 - x[i] +(1 + alpha*x[i])* rnorm(1)
        y[i, 2] <- 1 - x[i] - 0.5 * y[i, 1] + (1+alpha*x[i])*rnorm(1)
    }
}

xdim <- 3
X <- matrix(0, n, xdim)
X[,1] <- 1
X[,2] <- x
x2 <- rnorm(n)
X[,3] <- x2

## m <- 3
## n <- 10
## s <- rep(1:n,rep(m,n))
## x <- exp(rnorm(n*m))
## X <- cbind(1,x)
## u <- x*rnorm(m*n) + (1-x)*rf(m*n,3,3)
## a <- rep(rnorm(n),rep(m,n))
## y <- a + u
y1 <- y[, 1]

mod <- glm(R ~ y1 + X[, -1], family = "binomial")

eta <- fitted(mod)

pidi <- R * eta + (1 - R) * (1 - eta)

ystar <- y/pidi
Xstar <- X
Xstar[, -1] <- X[, -1]/pidi
ind <- c(seq(1, n * 2 - 1, by = 2), seq(2, n * 2, by = 2))

Xstarc <- rbind(Xstar, Xstar)
Xstarc2[ind, ] <- Xstarc

ystarc <- c(t(ystar))
Ry1 <- rep(1, n)
Rc <- c(rbind(Ry1, R))

library(quantreg)
mod <- rq(ystarc ~ Xstarc2[, -1], subset = Rc == 1)
coef <- coef(mod)
